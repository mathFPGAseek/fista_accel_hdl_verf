%--------------------------------------------------------------------------
% file: verify_1d_fft_float.m
% engr: rbd
% date : 11/3/24
% descr: Verify matching output of FPGA HDL after 1d FFT
%--------------------------------------------------------------------------
close all
clear all
clc           
%--------------------------------------------------------------------------
%% Run Vivado FPGA simulator; External tool
%--------------------------------------------------------------------------
disp(' Vivado 1D FFT should have been run ');
%--------------------------------------------------------------------------
%% Check 1-D FFT simulation
%--------------------------------------------------------------------------
convert_vectors_to_single; % after COPYING from viv wk "fft_1d_mem_raw_vevtors"

% Generate Xilinx MAT data to be checked
if exist('complex_image_array','var')
    disp(' Generate MAT Xilinx file for 1-D FFT checking');
    save(['./data/fft_float_' ...
        '1d_seq_matrix_fr_viv_sim.mat'],'complex_image_array');
    clearvars;
else
    error('Error: Could not generate MAT Xilinx file for testing');
end


% Generate Matlab model data to be checked against
run_xfft_v9_1_float_mex;

if exist('ImgByRow','var')
    disp(' Generate MAT file for 1-D FFT checking');
    save('./data/fft_float_1d_seq_matrix_fr_matlab.mat','ImgByRow'); %% no fftshift
    close all; clearvars;
else
    error('Error: Could not generate MAT file for testing');
end

% Perform 1-D FFT check
compare_hdl_sim_2_mat_sim_spec_1d_fft_float; % Checked Data & fftshift (Also reorders Viv vectors)

debug = 1;
