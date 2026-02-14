%-------------------------------------------------------------------------
%
% filename: convert_vectors_to_single.m
% engr: rbd
% date: 10/27/24
%
% descr:
% - raison d'etre:
% - Read in a file that is outputted from hardware simulation
% - so that we can compare against Matlab model
%
% - algo example:
% - temp3 = fi(-.0039,1,34,33);
% - temp3_bin = temp3.bin;
% - we get : '1111111110000000001101000110110111'
% - split into a fixed point string :
% - S = strsplit( '0.111111110000000001101000110110111', '.');
% - decimal positive if most sig = 0
%--------------------------------------------------------------------------
%% declarations
S_dot = ".";
S_array = {};
j = 0; % to make cell array  at zero
rows = 256; % equal column for square image
% - negative if most sig = 1
%fixed_point_length = 34;
float_point_length = 32;

% For debug mismatch use lines below, otherwise comment our and use
% fft1dmemevectors below these two lines
%rows = 32;
%fft1dmemvectors = importfile("C:\design\mig_zcu_104_2\fista_accel_verf\fista_accel_hdl_verf\fft_float\xfft_v9_1_bitacc_cmodel_nt64\debug\fft_1d_mem_raw_float_vectors_dbg_mismatch.txt", [1, Inf]); % TEMP DEBUG !!!


% Import Table
% DO NOT USE THIS LINE fft1dmemvectors = importfile("C:\design\mig_zcu_104_2\ddr4_0_ex\ddr4_0_ex.sim\sim_1\behav\xsim\fft_1d_mem_raw_vectors.txt", [1, Inf]);
%fft1dmemvectors  = importfile("C:\design\mig_zcu_104_2\ddr4_0_ex\ddr4_0_ex.sim\sim_1\behav\xsim\fft_1d_mem_raw_float_vectors.txt", [1, Inf]); % TEMP DEBUG !!!
% Added gen proc H logic
%fft1dmemvectors  = importfile("C:\design\mig_zcu_104_2022dot2\ddr4_0_ex_20250308_dev_work_with_v2022doot2.xpr\ddr4_0_ex\ddr4_0_ex.sim\sim_1\behav\xsim\fft_1d_mem_raw_float_vectors.txt", [1, Inf]); % TEMP DEBUG !!!
% Lenless Accel project cleaned up; regression test for 1d fft
% fft1dmemvectors  = importfile("C:\design\mig_zcu_104_2022dot2_clean\lensless_accel\lensless_accel.sim\sim_1\behav\xsim\fft_1d_mem_raw_float_vectors.txt",[1,Inf]);
% restart proj simulation 1/22/26 ; When we checked in 2/26 not matching
%fft1dmemvectors  = importfile("C:\design\mig_zcu_104_2022dot2_clean_recreate3\lensless_accel\lensless_accel.sim\sim_1\behav\xsim\fft_1d_mem_raw_float_vectors.txt",[1,Inf]);

% Capture fft1d @ memory tranpose module
fft1dmemvectors  = importfile("C:\design\mig_zcu_104_2022dot2_clean_recreate3\lensless_accel\lensless_accel.sim\sim_1\behav\xsim\MEM_TRANSPOSE_row_wr_mem_raw_fft_1d_float_vectors.txt",[1,Inf]);

% debug from 1d FFT with point source; 2nd and 3rd This has a bug with the
% StringSplit Array!!
%fft1dmemvectors = ["00001111011100011011000000010110100000111100001100000101010001110000";
%                   "00101111100100001000110111001011010010111101001011100111001101011100";];

size_string = size(fft1dmemvectors,1);



% Split into Real and Imag Strings
%[StringSplitArrayImagOut, StringSplitArrayRealOut] = SplitStringArray(fft1dmemvectors,fixed_point_length,size_string);
[StringSplitArrayImagOut, StringSplitArrayRealOut] = SplitStringArray(fft1dmemvectors,float_point_length,size_string); % TEMP DEBUG!!!

for i = 1: size_string
    sample_real_str = StringSplitArrayRealOut(i);
    sample_imag_str = StringSplitArrayImagOut(i);
    % debug
    %sample_str = "00010001000100011000100010001000"
    sample_real_hex = binToHex(sample_real_str); % char
    sample_imag_hex = binToHex(sample_imag_str); % char
    % convert binary to hex
    ArrayImageRealHex{i} = sample_real_hex; % store chars as cells to accumulate 
    ArrayImageImagHex{i} = sample_imag_hex; % store chars as cells to accumulate 

end


%?? convert hex to single
% This should give us our test vector ; Try this to make sure this works!
debug = 1;
%hexStr = '3F800000';
for i = 1 : size_string
    hexRealStr = ArrayImageRealHex{i};
    hexImagStr = ArrayImageImagHex{i};

    %debug 
    %hexRealStr  = '003F800000'; % 1 (dec)
    %hexRealStr  = '0015779688'; % 5E^-26
    floatRealValue(i) = hexToSingle(hexRealStr);
    floatImagValue(i) = hexToSingle(hexImagStr);
    debug = 1;
end 

debug = 1;


% convert to decimal - Pass String Array 
%[S_imag_array] = ConvStringArray2CellArray(StringSplitArrayImagOut,size_string);
%[S_real_array] = ConvStringArray2CellArray(StringSplitArrayRealOut,size_string);

% Write into an array to compare to 2-d fft model
r = size_string/rows;
c = r;
k = 1; % index for cell array
s_imag_numeric_array = zeros(r,c);
s_real_numeric_array = zeros(r,c);
for i = 1 : r
    for j = 1 : c
       %s_imag_numeric_array(i,j) = S_imag_array{k};
       s_imag_numeric_array(i,j) = floatImagValue(k);

       %s_real_numeric_array(i,j) = S_real_array{k};
       s_real_numeric_array(i,j) = floatRealValue(k);

       k = k +1;
    end
end

% Assign complex wt to imag components
s_imag_weighted_array = s_imag_numeric_array*1i;

% Create final complex array
complex_image_array = s_real_numeric_array + s_imag_weighted_array;
clearvars -except complex_image_array

debug = 1;
%--------------------------------------------------------------------
%% Functions
%--------------------------------------------------------------------
function hexStr = binToHex(binStr)

   %??? tomorrow try a nonzero hex
   %?? tried zero and got zer
   %?? why do I get a 1x2 string
   %?? I think that the binary string needs to be a char array

    binChar = convertStringsToChars(binStr);
    debug = 1;

    % Converts a binary string to hexadecimal
    %binStr - Input binary string (e.g., '110101101')
    %binStr = '110101101';
    % Ensure the length of the binary string is a multiple of 4
    % Padding with leading zeros if necessary
    n = length(binStr);
    if mod(n, 4) ~= 0
        binStr = [repmat('0', 1, 4 - mod(n, 4)), binStr];
    end
    
    % Convert binary string to decimal
    decimalValue = bin2dec(binStr(2));
    
    % Convert decimal to hexadecimal
    hexStr = dec2hex(decimalValue);
    
    % Display the result
    %disp(['Binary: ', binStr]);
    %disp(['Hexadecimal: ', hexStr]);
end

function floatValue = hexToSingle(hexStr)
    % Converts a hexadecimal string to a single-precision floating point number
    % hexStr - Hexadecimal string (e.g., '3F800000' for 1.0 in IEEE 754 format)
    
    % Convert hex string to a 32-bit unsigned integer
    uint32Value = uint32(hex2dec(hexStr));
    
    % Convert the 32-bit unsigned integer to single-precision floating point
    floatValue = typecast(uint32Value, 'single');
    
    % Display the result
    %disp(['Hexadecimal: ', hexStr]);
    %disp(['Single precision floating-point value: ', num2str(floatValue)]);
end

function [DecimalOut] = Conv2Dec(S)
    signV      = S{1} - '0';
    expV       = S{2} - '0';
    fracV      = S{3} - '0'; 
    signValue  = signV;
    expValue   = expV  * (2 .^ (numel(expV)-1:-1:0).'); 
    fracValue  = fracV * (2 .^ -(1:numel(fracV)).');
    if signValue == 0 % pos case
        DecimalOut = 2.^(expValue-127)*(1+fracValue);
    else
        DecimalOut = -1*2.^(expValue-127)*(1+fracValue);
    end
end
function fft1dmemvectors1_r1 = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  FFT1DMEMVECTORS1_R1 = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a string
%  array.
%
%  FFT1DMEMVECTORS1_R1 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  fft1dmemvectors1_r1 = importfile("C:\design\phd_ee\summer_2022_fpga_accel_admm_lenless_camera\admm_fpga_xilinx_project\ip\fft_256_test_2\xfft_0\cmodel\xfft_v9_1_bitacc_cmodel_nt64\fft_1d_mem_vectors.txt", [1, Inf]);
%
%  See also READMATRIX.
%
% Auto-generated by MATLAB on 10-Dec-2022 11:19:25

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "VarName1";
opts.VariableTypes = "string";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName1", "EmptyFieldRule", "auto");

% Import the data
fft1dmemvectors1_r1 = readmatrix(filename, opts);

end
function[S_array] = ConvStringArray2CellArray(fft1dmemvectors,size_string)
    % declarations
    S_dot = ".";
    S_array = {};
    j = 0; % to make cell array  at zero
    % convert to decimal
    for i = 1 : size_string 
    %for i = 32000 : size_string % TEMP DEBUG!!!
        S_string = fft1dmemvectors(i);
        S_dec = insertAfter(S_string,1,S_dot); % To separate sign from exp.
        S_dec = insertAfter(S_dec,10,S_dot); % To separate exp. from frac.
        S_format = strsplit(S_dec,'.');
        [DecimalOut] = Conv2Dec(S_format);
        S_array{j+1} = DecimalOut;
        j = j+1; 
    end
end
function[StringSplitArrayImagOut, StringSplitArrayRealOut] = SplitStringArray(StringIn,FixedPointLengthIn,SizeStringIn)
    Sdot = ".";
    k = 0;
    % create cell array
    CellArrayImag = {};
    CellArrayReal = {};
    % create string array
    StringSplitArrayImagOut = strings(SizeStringIn,1);
    StringSplitArrayRealOut = strings(SizeStringIn,1);

    for i = 1 : SizeStringIn
        Stemp = insertAfter(StringIn(i),FixedPointLengthIn,Sdot);
        Ssplit = strsplit(Stemp,Sdot);
        StringSplitImagTemp = Ssplit{1};
        StringSplitRealTemp = Ssplit{2};
        CellArrayImag{k+1} = StringSplitImagTemp;
        CellArrayReal{k+1} = StringSplitRealTemp;
        k = k+1;
    end

    % Create Character array
    CharSplitArrayImagOutTemp = char(CellArrayImag{1:SizeStringIn});
    CharSplitArrayRealOutTemp = char(CellArrayReal{1:SizeStringIn});

    % Create String array
    for i = 1 : SizeStringIn
        % Imag
        CharLineImagTemp = CharSplitArrayImagOutTemp(i,1:FixedPointLengthIn);
        StringLineImagTemp = convertCharsToStrings(CharLineImagTemp);
        StringSplitArrayImagOut(i) = StringLineImagTemp;
        %Real
        CharLineRealTemp = CharSplitArrayRealOutTemp(i,1:FixedPointLengthIn);
        StringLineRealTemp = convertCharsToStrings(CharLineRealTemp);
        StringSplitArrayRealOut(i) = StringLineRealTemp;
    end
end



