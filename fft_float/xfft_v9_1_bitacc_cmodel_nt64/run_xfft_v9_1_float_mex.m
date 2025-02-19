% % $RCSfile: run_xfft_v9_1_mex.m,v $ $Version: $ $Date: 2010/09/08 12:33:21 $
%
%  (c) Copyright 2008-2009 Xilinx, Inc. All rights reserved.
%
%  This file contains confidential and proprietary information
%  of Xilinx, Inc. and is protected under U.S. and
%  international copyright and other intellectual property
%  laws.
%
%  DISCLAIMER
%  This disclaimer is not a license and does not grant any
%  rights to the materials distributed herewith. Except as
%  otherwise provided in a valid license issued to you by
%  Xilinx, and to the maximum extent permitted by applicable
%  law: (1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND
%  WITH ALL FAULTS, AND XILINX HEREBY DISCLAIMS ALL WARRANTIES
%  AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, INCLUDING
%  BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-
%  INFRINGEMENT, OR FITNESS FOR ANY PARTICULAR PURPOSE; and
%  (2) Xilinx shall not be liable (whether in contract or tort,
%  including negligence, or under any other theory of
%  liability) for any loss or damage of any kind or nature
%  related to, arising under or in connection with these
%  materials, including for any direct, or any indirect,
%  special, incidental, or consequential loss or damage
%  (including loss of data, profits, goodwill, or any type of
%  loss or damage suffered as a result of any action brought
%  by a third party) even if such damage or loss was
%  reasonably foreseeable or Xilinx had been advised of the
%  possibility of the same.
%
%  CRITICAL APPLICATIONS
%  Xilinx products are not designed or intended to be fail-
%  safe, or for use in any application requiring fail-safe
%  performance, such as life-support or safety devices or
%  systems, Class III medical devices, nuclear facilities,
%  applications related to the deployment of airbags, or any
%  other applications that could lead to death, personal
%  injury, or severe property or environmental damage
%  (individually and collectively, "Critical
%  Applications"). Customer assumes the sole risk and
%  liability of any use of Xilinx products in Critical
%  Applications, subject only to applicable laws and
%  regulations governing limitations on product liability.
%
%  THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS
%  PART OF THIS FILE AT ALL TIMES. 
%-------------------------------------------------------------------
%
% Example code for FFT v8.0 MEX function
%
%-------------------------------------------------------------------
clf,clear all, close all
% Image to be tested
sizeX = 256;
sizeY = sizeX;
% Setup
ImgByRow = zeros(sizeX,sizeY);
PreImage = zeros(sizeX,sizeY);

% Generics for this smoke test
%generics.C_NFFT_MAX = 10;
generics.C_NFFT_MAX = 8;
%generics.C_NFFT_MAX = 256;

generics.C_ARCH = 1;
generics.C_HAS_NFFT = 0;
generics.C_USE_FLT_PT = 1;
generics.C_INPUT_WIDTH = 32; % Must be 32 if C_USE_FLT_PT = 1
generics.C_TWIDDLE_WIDTH = 24; % Must be 24 or 25 if C_USE_FLT_PT = 1
generics.C_HAS_SCALING = 0; % Set to 0 if C_USE_FLT_PT = 1
generics.C_HAS_BFP = 0; % Set to 0 if C_USE_FLT_PT = 1
generics.C_HAS_ROUNDING = 0; % Set to 0 if C_USE_FLT_PT = 1

channels = 1;

samples = 2^generics.C_NFFT_MAX;

% create a square full intensity at center of image

Max_value = 255; % uint8 rgb 
N_image = sizeX;
N_imp = 16;
[ image_out] = center_impulse(N_imp,N_image,Max_value);
%X = uint8(image_out);

% debug
load('C:\design\mig_zcu_104_2\fista_accel_verf\fista_accel_hdl_verf\fft_float\xfft_v9_1_bitacc_cmodel_nt64\data\point_src.mat'); 
%load('C:\design\mig_zcu_104_2\fista_accel_verf\fista_accel_hdl_verf\fft_float\xfft_v9_1_bitacc_cmodel_nt64\debug\hacked_img.mat'); % this gives us 1E-6
%load('C:\design\mig_zcu_104_2\fista_accel_verf\fista_accel_hdl_verf\fft_float\xfft_v9_1_bitacc_cmodel_nt64\debug\hacked_img_with_one_more_col_of_zeros.mat');
%load('C:\design\mig_zcu_104_2\fista_accel_verf\fista_accel_hdl_verf\fft_float\xfft_v9_1_bitacc_cmodel_nt64\debug\hacked_img_corrected.mat');

%X = hacked_img; % This gave us 1E-6
X = mem_array;

% Handle multichannel FFTs if required
for channel = 1:channels
   for i = 1 : sizeX

        % This is from point source using correct float reprsentation
        %input_raw = vpa(X(i,:));
         input_raw = X(i,:);
         input_raw = double(input_raw);

        %input_raw = double(double(X(i,:))/(Max_value+1)); % This gives us
        %1E-6

        % debug test vector
        %input_raw_sample = single(5.00000009770741e-26 + 5.00000009770741e-26i);
        %input_raw = repmat(input_raw_sample, [ 1 sizeX]);
        %input_raw = double(input_raw);


        % Create input data frame: constant data
        %constant_input = 0.5 + 0.5j;
        %input_raw(1:samples) = constant_input;
    
        if generics.C_USE_FLT_PT == 0
            % Set up quantizer for correct twos's complement, fixed-point format: one sign bit, C_INPUT_WIDTH-1 fractional bits
            % Alas, this function is not available in Octave, but the test values of 0.5+0.5j do not cause quantization errors anyway.
            %    q = quantizer([generics.C_INPUT_WIDTH, generics.C_INPUT_WIDTH-1], 'fixed', 'convergent', 'saturate');
            % Format data for fixed-point input
            %    input = quantize(q,input_raw);
            input = input_raw;
        else
            % Floating point interface - use data directly
            input = input_raw;
        end
  
        % Set point size for this transform
        %nfft = generics.C_NFFT_MAX;
        nfft = 256;
  
        % Set up scaling schedule: scaling_sch[1] is the scaling for the first stage
        % Scaling schedule to 1/N: 
        %    2 in each stage for Radix-4/Pipelined, Streaming I/O
        %    1 in each stage for Radix-2/Radix-2 Lite
        if generics.C_ARCH == 1 || generics.C_ARCH == 3
            scaling_sch = ones(1,floor(nfft/2)) * 2;
            if mod(nfft,2) == 1
                scaling_sch = [scaling_sch 1];
            end
        else
            scaling_sch = ones(1,nfft);
        end

        % Set FFT (1) or IFFT (0)
        %direction = 1;
        direction = 0;
      
        if channels > 1
            fprintf('Running the MEX function for channel %d...\n',channel)
        else
            fprintf('Running the MEX function...\n')      
        end

        % debug
        if i == 121
            debug = 1;
            %input = [input(1:120),input(120),input(121:128),input(129:end-1) ];
            %nfft = 256;
        end
  
        % Run the MEX function
        [output, blkexp, overflow] = xfft_v9_1_bitacc_mex(generics, nfft, input, scaling_sch, direction);
  
        % Gather up 1D FFT for image
        %ImgByRow(i,:) = (fftshift(output));
         ImgByRow(i,:) = output;


   end % input lines

  debug = 1;

  % Check output values are correct
  % The FFT of constant input data is an impulse
  % Therefore all output samples should be zero except for the first
  % The value of the first sample depens on the type of scaling used
  %{
  if generics.C_USE_FLT_PT == 0
      if generics.C_HAS_SCALING == 0
        expected_xk_re_0 = constant_input * (2^generics.C_NFFT_MAX);
      else
        expected_xk_re_0 = constant_input;
      end
  else 
     expected_xk_re_0 = 512 + 512j; 
  end

  % Check xk_re and xk_im data: Only xk_re[0] should be non-zero  
  if output(1) ~= expected_xk_re_0
      if channels > 1
        error('ERROR: Channel %d xk_re[0] is incorrect: expected %f + j%f, actual %f + j%f\n',channel,real(expected_xk_re_0),imag(expected_xk_re_0),real(output(1)),imag(output(1)))
      else
        error('ERROR: xk_re[0] is incorrect: expected %f + j%f, actual %f + j%f\n',real(expected_xk_re_0),imag(expected_xk_re_0),real(output(1)),imag(output(1)))          
      end
  end
  
  % Check all other sample values are zero
  for n = 2:samples
    if output(n) ~= 0 + 0j
        if channel > 1
          error('ERROR: Channel %d output sample %d is incorrect: expected %f +j%f, actual %f + j%f\n',channel,n,0.0,0.0,real(output(1)),imag(output(1)))
        else 
          error('ERROR: output sample %d is incorrect: expected %f +j%f, actual %f + j%f\n',n,0.0,0.0,real(output(1)),imag(output(1)))
        end
    end
  end
  %}
  % Check if blkexp used: should be nfft
  if generics.C_HAS_BFP == 1
    if blkexp ~= nfft
        if channels > 1
          error('ERROR: Channel %d blkexp is incorrect.  Expected value %d\n',channel,nfft)
        else
          error('ERROR: blkexp is incorrect.  Expected value %d\n',nfft)            
        end
    end
  end
  
  % Check overflow if used: scaling schedule should ensure that overflow never occurs
  if generics.C_HAS_SCALING == 1 && generics.C_HAS_BFP == 0
    if overflow == 1
        if channels > 1
          error('ERROR: Channel %d overflow is incorrect\n',channel)
        else
          error('ERROR: overflow is incorrect\n')            
        end
    end
  end
  
  % If we got to here, simulation outputs are correct
  if channels > 1
      fprintf('Test completed successfully\n',channel)
  else
      fprintf('Test completed successfully\n')      
  end
  
end

%--------------------------------------------------
%% functions
%--------------------------------------------------


function [ image_out] = center_impulse(N_imp,N_image,Max_value)
% center a N_impxN_imp inside of a N_image x N_image
% Only use square images that are 2^N( for ease use of FFT)
%X1 = zeros(N_image,N_image);
X1(1:N_image,1:N_image) = 5.8775e-39;
delta_row_lower = (N_image - N_imp)/2;
row_upper_start_pixel = N_image - delta_row_lower - N_imp;
row_upper_end_pixel   = row_upper_start_pixel + N_imp;
row_lower_start_pixel = N_image - delta_row_lower;
row_lower_end_pixel   = row_lower_start_pixel + N_imp;

% default to rgb uint8
for i = row_upper_start_pixel : row_upper_end_pixel
   for j = row_upper_start_pixel : row_upper_end_pixel
       X1(i,j) = Max_value;
   end
end

debug = 1;
% debug
lims = [ -256 256];
clim = [ 0 256];
figure(1), clf
plot(1)
imageh = imagesc(zeros(N_image));
axis square, axis off, axis xy
set(gca,'xlim',[lims(2)-30 lims(2)+30],'ylim',[lims(2)-30 lims(2)+30],'clim',[clim(1) clim(2)])
title('Test Centered Impulse')

set(imageh,'CData',X1);

debug = 1;
image_out = X1;


end



