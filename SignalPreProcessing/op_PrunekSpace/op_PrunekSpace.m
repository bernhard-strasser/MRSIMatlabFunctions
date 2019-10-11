function [MRStruct] = op_PrunekSpace(MRStruct, Settings)
%
% op_PermuteMRData Prune k-Space of MR Data to Size Given in Header
%
% This function was written by Bernhard Strasser, Oct 2019.
%
%
% The function removes the zeros that are sometimes present in (MRSI) raw data. These zeroes are there because the Line (kx) and Phase (ky) indices in the mdh are on
% purpose wrong. E.g. if you have a matrix size of 12x12, the k-space center will have Line and Phase indices of 9x9 instead of 7x7. This is in preparation because
% the (MRSI) data is zerofilled to the next power of 2 (16x16 in this example) in ICE. However, usually we don't like this zerofilling in our reco. 
% This function gets rid of that extra zeros.
%
%
% [MRStruct] = op_PermuteMRData(MRStruct,Settings)
%
% Input: 
% -         MRStruct               ...      The structure containing the MR data. Can have different fields. For this function, it
%                                           it has to have field
%                                           * .Data
% -         Settings               ...      Structure with settings how data should be processed. Should have field:
%                                           * .PermuteVec: Data will be permuted to according to this vector.
%
% Output:
% -         MRStruct               ...      The structure containing the MR data. Can have different fields.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 



%% Prep
if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isfield(Settings,'RemoveDataFrom'))
    Settings.RemoveDataFrom = 'Beginning';
end
if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct.Par,'DataSize'))
    MRStruct.Par.DataSize = size(MRStruct.Data);
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end


%% Prune Data

TargSize = [MRStruct.RecoPar.nPhasEnc MRStruct.RecoPar.nFreqEnc MRStruct.RecoPar.nPartEnc*MRStruct.RecoPar.nSLC];

if(strcmp(Settings.RemoveDataFrom,'Beginning'))
    StartInd = MRStruct.RecoPar.DataSize(1:3) - TargSize + 1;
end
MRStruct.Data = MRStruct.Data(StartInd(1):end,StartInd(2):end,StartInd(3):end,:,:,:,:,:,:,:,:,:,:,:);

if(isfield(MRStruct,'NoiseData'))
    MRStruct.NoiseData = MRStruct.NoiseData(StartInd(1):end,StartInd(2):end,StartInd(3):end,:,:,:,:,:,:,:,:,:,:,:);
end
MRStruct.RecoPar.DataSize = size(MRStruct.Data);


%% Postparations
MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



