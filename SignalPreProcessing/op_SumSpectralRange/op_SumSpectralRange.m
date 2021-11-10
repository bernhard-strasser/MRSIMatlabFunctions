function [MRStructOut,Settings] = op_SumSpectralRange(MRStruct, Settings)
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
%                                           * .SumSpectralRange_ppm: Range for summing.
%                                           * .TakeRealAbsImagComplex: Function to call before summing, .e.g @abs, @real, @imag etc.
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

%% Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'SumSpectralRange_ppm'))
    Settings.SumSpectralRange_ppm{1} = [-999 999];  % Everything
end
if(~iscell(Settings.SumSpectralRange_ppm))
   Tmp = Settings.SumSpectralRange_ppm;
   Settings = rmfield(Settings,'SumSpectralRange_ppm');
   Settings.SumSpectralRange_ppm{1} = Tmp;
end

if(~isfield(Settings,'TakeRealAbsImagComplex'))
    Settings.TakeRealAbsImagComplex = @abs;
end


if(isfield(MRStruct,'RecoPar'))
    CurPar = MRStruct.RecoPar;
else
    CurPar = MRStruct.Par;    
end




%% Prepare Data

chemy = compute_chemshift_vector(CurPar);
MRStruct.Data = feval(Settings.TakeRealAbsImagComplex,fftshift(fft(MRStruct.Data,[],4),4));
MRStructOut = MRStruct; MRStructOut.Data = [];



%% Sum Data

for ii = 1:numel(Settings.SumSpectralRange_ppm)
    
    f_first = FindClosestIndex(chemy,max(Settings.SumSpectralRange_ppm{ii})); Settings.f_first{ii} = f_first{1};
    f_end = FindClosestIndex(chemy,min(Settings.SumSpectralRange_ppm{ii})); Settings.f_end{ii} = f_end{1};
    
    MRStructOut.Data = cat(4,MRStructOut.Data,sum(MRStruct.Data(:,:,:,Settings.f_first{ii}:Settings.f_end{ii},:,:,:,:,:),4));
end



%% Postparations

MRStructOut = supp_UpdateRecoSteps(MRStructOut,Settings);



