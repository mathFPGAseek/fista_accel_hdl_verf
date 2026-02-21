%% fista_deblur_corner_conventions_demo.m
% Goal:
%   Remove “fear” by showing a complete deblur + partial measurement example where:
%     - PSF and measurement are created/displayed in a *centered* (fftshift) way
%     - BUT the actual FISTA operators A and A* are *FPGA-friendly*:
%         * no fftshift/ifftshift inside iterations
%         * “corner” indexing (DC at (1,1))
%         * corner-crop measurement (top-left My×Mx)
%         * corner-pad in A*
%   At the end, we only fftshift for DISPLAY (“WAH LA”).

clear; clc; close all;
rng(0);

%% ------------------------ Sizes ------------------------
Ny = 32; Nx = 32;        % unknown image size (FFT size in FPGA)
My = 20; Mx = 20;        % measurement window (corner crop)

assert(My<=Ny && Mx<=Nx);

%% ------------------------ Build a ground-truth image (for demo) ------------------------
% Make something with structure so deblur looks meaningful.
x0_centered = zeros(Ny,Nx);
x0_centered( 9:12,  8:24) = 1.0;     % bar
x0_centered(18:24, 10:13) = 0.7;     % vertical feature
x0_centered(22:26, 20:27) = 0.9;     % block
% Add a “point”
x0_centered(16,16) = 1.5;

% This x0_centered is in a "display centered" coordinate sense.
% Convert ONCE into the FPGA/corner layout for the actual optimization variable:
x0 = ifftshift(x0_centered);

%% ------------------------ Build a centered PSF (for demo) ------------------------
% Centered Gaussian PSF (peak in the middle for visualization)
sigma = 2.0;
[xx,yy] = meshgrid( (-floor(Nx/2)):(ceil(Nx/2)-1), (-floor(Ny/2)):(ceil(Ny/2)-1) );
psf_centered = exp(-(xx.^2 + yy.^2)/(2*sigma^2));
psf_centered = psf_centered / sum(psf_centered(:));  % normalize energy

% Convert ONCE to an uncentered OTF H suitable for corner FFT usage:
% IMPORTANT: this is the correct “physics” conversion: H = FFT{psf} with PSF at (1,1)
H = fft2(ifftshift(psf_centered));   % uncentered frequency response (FPGA-friendly)

%% ------------------------ Define corner-crop/pad ------------------------
corner_crop = @(z) z(1:My, 1:Mx);
corner_pad  = @(w) pad_to_corner(w, Ny, Nx);

%% ------------------------ Define forward A and adjoint A* (NO SHIFTS) ------------------------
% Forward blur (circular conv via FFT) then corner-crop measurement
A = @(x) corner_crop( ifft2( H .* fft2(x) ) );

% Adjoint: corner-pad back to full grid, then apply conjugate frequency response
Aadj = @(y) ifft2( conj(H) .* fft2( corner_pad(y) ) );

%% ------------------------ Simulate a centered measurement (like many MATLAB demos) ------------------------
% Create “full blurred image” on the full grid (in corner layout)
blur_full = ifft2( H .* fft2(x0) );    % no shifts

% Measurement is corner-cropped (FPGA measurement convention)
b = corner_crop(blur_full);

% Add a little noise
noise_sigma = 0.01;
b = b + noise_sigma*(randn(My,Mx) + 1j*0);  % keep real-ish for easy viewing

% For DISPLAY ONLY: show a centered version of the measurement
% (If you simulated a “centered measurement” in some other code, this is the place
%  you’d typically do fftshift for display, NOT for the algorithm.)
b_centered_display = fftshift(b);  % just to show something "center-ish" on plots

%% ------------------------ FISTA setup (L2 data + L1 prior) ------------------------
% Solve: minimize 0.5||A x - b||_2^2 + lambda * ||x||_1
% (L1 is just to make it a proper “prox” demo; you can set lambda=0 for plain LS.)

lambda = 0.01;         % adjust to taste
nIter  = 200;

% Step size t = 1/L, where L >= ||A||^2
% For convolution with H and a crop, a safe bound is:
L = max(abs(H(:)).^2);     % crop does not increase operator norm
t = 1 / L;

% FISTA variables
x  = zeros(Ny,Nx);         % start at 0 (corner layout)
z  = x;
tk = 1;

% Soft threshold for complex values (magnitude shrink)
soft = @(u,th) soft_complex(u, th);

%% ------------------------ Run FISTA (NO SHIFTS ANYWHERE) ------------------------
for k = 1:nIter
    r = A(z) - b;                 % residual in measurement domain (My×Mx)
    g = Aadj(r);                  % gradient in image domain (Ny×Nx)

    x_new = soft(z - t*g, lambda*t);

    tk_new = (1 + sqrt(1 + 4*tk^2))/2;
    z = x_new + ((tk - 1)/tk_new)*(x_new - x);

    x = x_new;
    tk = tk_new;

    % Optional: simple progress monitor
    if mod(k,50)==0
        f = 0.5*norm(r(:))^2 + lambda*sum(abs(x(:)));
        fprintf('iter %3d: objective %.6e\n', k, f);
    end
end

x_hat = x;  % reconstructed in CORNER layout

%% ------------------------ “WAH LA”: shift only for display ------------------------
x_hat_centered_display = fftshift(x_hat);

%% ------------------------ Plots ------------------------
figure('Name','Centered display (human-friendly)'); colormap gray;

subplot(2,3,1);
imagesc(x0_centered); axis image off; title('True x (centered display)'); colorbar;

subplot(2,3,2);
imagesc(psf_centered); axis image off; title('PSF (centered display)'); colorbar;

subplot(2,3,3);
imagesc(real(b_centered_display)); axis image off; title('Measured b (display-only)'); colorbar;

subplot(2,3,4);
imagesc(abs(fftshift(blur_full))); axis image off; title('Full blurred |Ax| (display-only)'); colorbar;

subplot(2,3,5);
imagesc(real(x_hat_centered_display)); axis image off; title('Reconstruction (centered display)'); colorbar;

subplot(2,3,6);
imagesc(real(x_hat_centered_display - x0_centered)); axis image off; title('Error (centered display)'); colorbar;

sgtitle('All ITERATIONS used corner FFT + corner crop/pad (FPGA style). fftshift used ONLY for plots.');

%% ------------------------ Adjoint sanity check (optional) ------------------------
ip = @(u,v) sum(conj(u(:)).*v(:));  % <u,v> = u^H v

xt = randn(Ny,Nx) + 1j*randn(Ny,Nx);
yt = randn(My,Mx) + 1j*randn(My,Mx);

lhs = ip(A(xt), yt);
rhs = ip(xt, Aadj(yt));
fprintf('\nAdjoint check: |<A x,y> - <x,A* y>| = %.3e\n', abs(lhs-rhs));

%% ------------------------ Local functions ------------------------
function X = pad_to_corner(Y, Ny, Nx)
    X = zeros(Ny, Nx, 'like', Y);
    [My,Mx] = size(Y);
    X(1:My, 1:Mx) = Y;
end

function y = soft_complex(x, th)
% Complex soft-thresholding: shrink magnitude, preserve phase.
    ax = abs(x);
    scale = max(ax - th, 0) ./ max(ax, eps(class(ax)));
    y = x .* scale;
end
