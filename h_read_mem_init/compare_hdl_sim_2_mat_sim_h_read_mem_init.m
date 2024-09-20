%--------------------------------------------------------------------------
% file: compare_hdl_sim_2_mat_sim_spec_1d_fft.m
% engr: rbd
% date : 12/21/22
% raison d'etre: match vectors hdl sim to matlab
% descr/instrs:
%--------------------------------------------------------------------------
%clc
%clear all
%close all

% rename arrays
ImgByColFrHdlSimSeq = complex_image_array_hdl; % import fft-2d ordered and shifted from hdl sim
ImgByColFrMatSimSeq = complex_image_array_mat; % import fft-2d  ordered and shifted from matlab

clear vars complex_image_array PreImage

re_HdlSim = real(ImgByColFrHdlSimSeq(:));
im_HdlSim = imag(ImgByColFrHdlSimSeq(:));

re_MatSim = real(ImgByColFrMatSimSeq(:));
im_MatSim = imag(ImgByColFrMatSimSeq(:));

figure(1)
%title('Real differences HDL Sim vs Matlab model');
plot(re_HdlSim - re_MatSim);
title('Real differences HDL Sim vs Matlab model');


figure(2)
%title('Imaginary differences HDL Sim vs Matlab model ');
plot(im_HdlSim - im_MatSim );
title('Imaginary differences HDL Sim vs Matlab model ');


disp('Testing complete, check plots  ')

debug = 1;
