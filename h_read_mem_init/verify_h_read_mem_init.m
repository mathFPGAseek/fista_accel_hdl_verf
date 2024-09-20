%--------------------------------------------------------------------------
% file: verify_h_read_mem_init.m
% engr: rbd
% date : 9/2/24
% descr: Verify matching output of FPGA HDL after H read Init
%
% Read in MAT file: fft_2d_A_fr_mat generated from
% ...verify_2d_fft:
% ...File : compare_hdl_sim_2_mat_sim_spec_2d_fft.m  
% ...signal: fft_2d_A_mem.reorderedspectrumMatSim
%
%--------------------------------------------------------------------------
close all
clear all
clc           

%--------------------------------------------------------------------------
%% Run Vivado FPGA simulator; External tool+
%--------------------------------------------------------------------------
disp(' Vivado 2D H read mem init should have been run ');
%--------------------------------------------------------------------------
%% Check H read mem init simulation
%--------------------------------------------------------------------------
convert_vectors_to_decimal; % after COPYING from viv wk

% Generate Xilinx MAT data to be checked
if exist('complex_image_array','var')
    disp(' H read mem init Exists for checking');  
else
    error('Error: Could not find H init mem file for testing');
end

% rename 
complex_image_array_hdl = complex_image_array;
clearvars -except complex_image_array_hdl;


% Generate Matlab model data to be checked against
load('./data/fft_2d_A_fr_mat.mat');

if exist('reorderedspectrumMatSim','var')
    disp(' Generate MAT file for H Init checking');
else
    error('Error: Could not find var from file for testing');
end

% rename 
complex_image_array_mat = reorderedspectrumMatSim;
clearvars -except complex_image_array_hdl complex_image_array_mat;


% Perform H read mem init check;
compare_hdl_sim_2_mat_sim_h_read_mem_init; 
debug = 1;
