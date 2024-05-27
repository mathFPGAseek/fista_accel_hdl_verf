%--------------------------------------------------------------------------
% file: verify_2d_fft.m
% engr: rbd
% date : 4/27/24
% descr: Verify matching output of FPGA HDL after 2d FFT
%--------------------------------------------------------------------------
close all
clear all
clc           

%--------------------------------------------------------------------------
%% Run Vivado FPGA simulator; External tool
%--------------------------------------------------------------------------
disp(' Vivado 2D FFT should have been run ');
%--------------------------------------------------------------------------
%% Check 2-D FFT simulation
%--------------------------------------------------------------------------
convert_vectors_to_decimal; % after COPYING from viv wk "fft_2d_mem_raw_vevtors"

% Generate Xilinx MAT data to be checked
if exist('complex_image_array','var')
    disp(' Generate MAT Xilinx file for 2-D FFT checking');
    save('./data/fft_2d_seq_matrix_fr_viv_sim.mat','complex_image_array');
    clearvars;
else
    error('Error: Could not generate MAT Xilinx file for testing');
end


% Generate Matlab model data to be checked against
run_xfft_v9_1_pipe_mex_2d_fft_no_shift;

if exist('PreImage','var')
    disp(' Generate MAT file for 2-D FFT checking');
    save('./data/fft_2d_seq_matrix_fr_matlab.mat','PreImage'); %% no fftshift on 2nd fft
    close all; clearvars;
else
    error('Error: Could not generate MAT file for testing');
end

% Perform 2-D FFT check; Note this is not adjusted for FFT scaling
% This is just a comparison between FFT in Vivado and FFT in Matlab
compare_hdl_sim_2_mat_sim_spec_2d_fft; % Checked Data & fftshift (Also reorders Viv vectors)

debug = 1;
