% $RCSfile: run_xfft_v9_1_mex.m,v $ $Version: $ $Date: 2010/09/08 12:33:21 $
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

% Generics for this smoke test matching FFT C model that we synthesized
generics.C_NFFT_MAX = log2(sizeX);
generics.C_ARCH = 3;
generics.C_HAS_NFFT = 0;
generics.C_USE_FLT_PT = 0;
generics.C_INPUT_WIDTH = 34; % Must be 32 if C_USE_FLT_PT = 1
generics.C_TWIDDLE_WIDTH = 34; % Must be 24 or 25 if C_USE_FLT_PT = 1
generics.C_HAS_SCALING = 1; % Set to 0 if C_USE_FLT_PT = 1
generics.C_HAS_BFP = 0; % Set to 0 if C_USE_FLT_PT = 1
generics.C_HAS_ROUNDING = 0; % Set to 0 if C_USE_FLT_PT = 1

channels = 1;

samples = 2^generics.C_NFFT_MAX;

% create a square full intensity at center of image

Max_value = 255; % uint8 rgb 
N_image = sizeX;
N_imp = 16;
[ image_out] = center_impulse(N_imp,N_image,Max_value);
X = uint8(image_out);

      
if channels > 1
     fprintf('Running the MEX function for channel %d...\n',channel)
else
     fprintf('Running the MEX function...\n')      
end

% Handle multichannel FFTs if required
for channel = 1:channels
 % Process FFT Rows first
 tic
 for i = 1 : sizeX
  input = double(double(X(i,:))/(Max_value+1)); 

  % Set point size for this transform
  nfft = generics.C_NFFT_MAX;

  % debug
  if ( i == 120)
      debug = 1;
  end
  
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
  direction = 0;      
  
  % Run the MEX function
  [output, blkexp, overflow] = xfft_v9_1_bitacc_mex(generics, nfft, input, scaling_sch, direction);
  ImgByRow(i,:) = (fftshift(output));
  %ImgByRow(i,:) = output;

  % Check overflow if used: scaling schedule should ensure that overflow never occurs
  if generics.C_HAS_SCALING == 1 && generics.C_HAS_BFP == 0
    if overflow == 1
        if channels > 1
          error('ERROR: Channel %d overflow, row %d is incorrect\n',channel,i)
        else
          error('ERROR: overflow, row %d is incorrect\n',i)            
        end
    end
  end

 end % for parsing rows
 toc
 tic
 % If we got to here, simulation outputs are correct
 if channels > 1
      fprintf('Test completed, Channel %d,1st FFT stage successfully\n',channel)
 else
      fprintf('Test completed 1st FFT stage successfully\n')      
 end
 % Process FFT by Col
 for j = 1 : sizeY
     XCol = ImgByRow(:,j);
     input = XCol;
     % Set point size for this transform
     nfft = generics.C_NFFT_MAX;
  
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
     direction = 1;      
  
     % Run the MEX function
     [output, blkexp, overflow] = xfft_v9_1_bitacc_mex(generics, nfft, input, scaling_sch, direction);
     %PreImage(:,j) = (fftshift(output));
      PreImage(:,j) = output;


  % Check overflow if used: scaling schedule should ensure that overflow never occurs
  if generics.C_HAS_SCALING == 1 && generics.C_HAS_BFP == 0
    if overflow == 1
        if channels > 1
          error('ERROR: Channel %d overflow, col %d is incorrect\n',channel,j)
        else
          error('ERROR: overflow, col %d is incorrect\n',j)            
        end
    end
  end

 end % for parsing cols
 toc
 % If we got to here, simulation outputs are correct
 if channels > 1
      fprintf('Test completed, Channel %d,2nd FFT stage successfully\n',channel)
 else
      fprintf('Test completed 2nd FFT stage successfully\n')      
 end
     
end % for channels

% Plot results
img = abs(PreImage.*2048.*256);  % multiply adjust 2048 for shifts thru 2 FFTs
                            % mulitple adjust by 256 for original
                            % normalization
figure(1), clf
plot(1)
imageh = imagesc(zeros(sizeX));

debug = 1;

% For displaying image intensity of ffts
Clim1 = 100; 
createfigure(img,Clim1)

debug = 1;

%--------------------------------------------------
%% functions
%--------------------------------------------------

function [ image_out] = center_impulse(N_imp,N_image,Max_value)
% center a N_impxN_imp inside of a N_image x N_image
% Only use square images that are 2^N( for ease use of FFT)
X1 = zeros(N_image,N_image);
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

function createfigure(cdata1,Clim1)
%CREATEFIGURE(cdata1)
%  CDATA1:  image cdata

%  Auto-generated by MATLAB on 03-Sep-2022 13:09:48

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
axis off
hold(axes1,'on');

% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');

% Create title
%title('Amplitude Spectrum w/impulse square size',num2str(N_imp));
title('Amplitude Spectrum w/impulse square size');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[-7.20655734386856 296.061659755709]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-23.8000217425869 279.46819535699]);
box(axes1,'on');
axis(axes1,'square');
%hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'CLim',[0 Clim1],'Layer','top');

hold(axes1,'off');
end



