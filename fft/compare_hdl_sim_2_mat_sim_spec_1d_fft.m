%--------------------------------------------------------------------------
% file: compare_hdl_sim_2_mat_sim_spec_1d_fft.m
% engr: rbd
% date : 12/21/22
% raison d'etre: match vectors hdl sim to matlab
% descr/instrs:
% Rearrange sims from hdl sim and matlab to compare.
% Steps :
% 0. Copy fft_1d_mem_raw_vectors.txt to this working directory.
% 1. Run convert_vectors_to_decimal.m ... This converts the HDL sim text
% file in binary to decimal values in an array called complex_image_array
% 2. Copy 'complex_image_array' to mat file called:
% 'fft_1d_seq_matrix_fr_viv_sim.mat'
% 3. Run run_xfft_v9_1_pipe_mex_1d_fft_no_shift.m up to generation of rows
% with no shift
% 4. Run clearvars -except ImgByRow
% 5. Copy ImgByRow to mat file called:
% 'fft_1d_seq_matrix_fr_matlab.mat'
% 6. Run this file
% Expected Results:
% HDL is sequenced from 1 to 256
% Matlab is sequenced from 1 to 256 because no fftshift
% error is less than 2 *10^-3
%--------------------------------------------------------------------------
clc
clear all
close all

% load hdl sim and matlab sim
load('./data/fft_1d_seq_matrix_fr_viv_sim.mat');
load('./data/fft_1d_seq_matrix_fr_matlab.mat');

% rename arrays
ImgByRowFrHdlSimSeq = complex_image_array; % import fft-1d from hdl sim
ImgByRowFrMatSimSeq = ImgByRow; % import fft-1d from matlab

clear vars complex_image_array ImgByRow

% reorder array conversion from hdl sim

hdl_sim_order = [0  128 64 192	32	160	96	224	16	144	80	208	48	176	112	240	8	136	72	200	40	168	104	232	24	152	88	216	56	184	120	248 ...
                 4	132 68 196	36	164	100 228	20	148	84	212	52	180	116	244	12	140	76	204	44	172	108	236	28	156	92	220	60	188	124	252 ...
                 2	130 66 194	34	162	98	226	18	146	82	210	50	178	114	242	10	138	74	202	42	170	106	234	26	154	90	218	58	186	122	250 ...
                 6	134 70 198	38	166	102 230	22	150	86	214	54	182	118	246	14	142	78	206	46	174	110	238	30	158	94	222	62	190	126	254 ...
                 1	129 65 193	33	161	97	225	17	145	81	209	49	177	113	241	9	137	73	201	41	169	105	233	25	153	89	217	57	185	121	249 ...
                 5  133 69 197	37	165	101 229	21	149	85	213	53	181	117	245	13	141	77	205	45	173	109	237	29	157	93	221	61	189	125	253 ...
                 3	131 67 195	35	163	99	227	19	147	83	211	51	179	115	243	11	139	75	203	43	171	107	235	27	155	91	219	59	187	123	251 ...
                 7	135 71 199	39	167	103 231	23	151	87	215	55	183	119	247	15	143	79	207	47	175	111	239	31	159	95	223	63	191	127	255];

% init
rows = size(ImgByRowFrHdlSimSeq,1);
col = rows;
reorderedHdlSim = zeros(rows,col);


%% reorder hdlsim as sequence
k = 1;
index = 0;
for i = 1 : rows
    for j = 1 : col
        index = hdl_sim_order(k);
        index = index + 1; % matlab does not like zero
        value_to_write = ImgByRowFrHdlSimSeq(i,j);
        reorderedHdlSim(i,index) = value_to_write;
        k = k + 1;
    end
    k = 1;
end

%% reorder as a spectrum that has center frquency as dc( bin 0)
reorderedHdlSim         = fftshift(reorderedHdlSim,2);
reorderedspectrumMatSim = fftshift(ImgByRowFrMatSimSeq,2);

%{
% arrange sequences as columns for debugging
TransposeImgByRowFrHdlSimSeq     = reorderedHdlSim(120,:)';
TransposeImgByRowFrMatSimSeq     = reorderedspectrumMatSim(120,:)';


% compare mag
reorderedMatSimMag  = abs(reorderedspectrumMatSim );
reorderedHdlSimMag  = abs(reorderedHdlSim);

diff = reorderedMatSimMag  - reorderedHdlSimMag;
surf(diff,'edgecolor','none');
%}
re_reorderedHdlSim = real(reorderedHdlSim(:));
im_reorderedHdlSim = imag(reorderedHdlSim(:));

re_reorderedspectrumMatSim = real(reorderedspectrumMatSim(:));
im_reorderedspectrumMatSim = imag(reorderedspectrumMatSim(:));

figure(1)
title('Real differences HDL Sim vs Matlab model');
plot(re_reorderedHdlSim - re_reorderedspectrumMatSim);

figure(2)
title('Imaginary differences HDL Sim vs Matlab model ');
plot(im_reorderedHdlSim - im_reorderedspectrumMatSim );

disp('Testing complete, check plots  ')

debug = 1;
