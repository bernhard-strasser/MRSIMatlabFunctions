function [MRStruct,AdditionalOut] = op_NoiseDecorrelateData(MRStruct, NoiseStructOrNoiseCorrMat, Settings)
%
% op_PermuteMRData Permute MR Data
%
% This function was written by Bernhard Strasser, Oct 2019.
%
%
% The function just runs the "permute" function on MR data as specified by Settings.PermuteVec.
% It also takes care of things like updating the .Par.DataSize, and permuting the NoiseData, if available.
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

if(~isstruct(MRStruct) || numel(fieldnames(MRStruct)) < 1)
    return;
end

if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isfield(Settings,'CreateNoiseData'))
    Settings.CreateNoiseData = true;
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


%% Calculate NoiseCorrMat


if(isstruct(NoiseStructOrNoiseCorrMat))
    
    Noise_mat = reshape(NoiseStructOrNoiseCorrMat.Data,[size(NoiseStructOrNoiseCorrMat.Data,1) numel(NoiseStructOrNoiseCorrMat.Data)/size(NoiseStructOrNoiseCorrMat.Data,1)]);
    Noise_mat = Noise_mat(:,20:end);
    NoiseScalingFactor = sqrt(NoiseStructOrNoiseCorrMat.Par.Dwelltimes(1) / MRStruct.Par.Dwelltimes(1));
    % NoiseScalingFactor = 1;
    Noise_mat = Noise_mat * NoiseScalingFactor;
    NoiseCorrMat = 1/(size(Noise_mat,2)) * (Noise_mat * Noise_mat');
else
    NoiseCorrMat = NoiseStructOrNoiseCorrMat;
end
AdditionalOut.NoiseCorrMat = NoiseCorrMat;


%% Uncell MRStruct data if necessary
    
WasCell_flag = false;
if(iscell(MRStruct.Data))
    MRStruct.Data = cat(1,MRStruct.Data{:});
    WasCell_flag = true;
end


%% Perform NoiseDecorrelation

MRStruct = op_PermuteMRData(MRStruct,[6 1 2 3 4 5 7:numel(size(MRStruct.Data))]);
MRStruct.Data = PerformNoiseDecorrelation(MRStruct.Data,NoiseCorrMat);
MRStruct = op_PermuteMRData(MRStruct,[2 3 4 5 6 1 7:numel(size(MRStruct.Data))]);


%% Cell MRStruct data if necessary

if(WasCell_flag)
    CurPt = 1;
    Bak = cell([1 numel(MRStruct.RecoPar.TrajPts)]);
    for CurAI = 1:numel(MRStruct.RecoPar.TrajPts)
        Bak{CurAI} = MRStruct.Data(CurPt:CurPt+MRStruct.RecoPar.TrajPts(CurAI)-1,:,:,:,:,:,:,:);
        CurPt = CurPt + MRStruct.RecoPar.TrajPts(CurAI);
    end
    MRStruct.Data = Bak; clear Bak;
end



%% 
if(Settings.CreateNoiseData)
    if(WasCell_flag)
        for CurAI = 1:numel(MRStruct.RecoPar.TrajPts)
            MRStruct.NoiseData{CurAI} = randn(size(MRStruct.Data{CurAI})); 
        end
    else
        MRStruct.NoiseData = randn(size(MRStruct.Data)); 
    end
end


%% Postparations
MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



