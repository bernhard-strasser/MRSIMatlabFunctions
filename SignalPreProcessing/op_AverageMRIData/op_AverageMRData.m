function [MRStruct,AdditionalOut] = op_AverageMRData(MRStruct, Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%%

if(~exist('Settings','var'))
   Settings = struct; 
    end
if(~isfield_recursive(Settings,'UndoWeightedAveraging_flag'))
    Settings.UndoWeightedAveraging_flag = true;
end
if(~isfield(MRStruct,'RecoPar'))
	MRStruct.RecoPar = MRStruct.Par;
    end

if(isfield(MRStruct.RecoPar,'dimnames'))
    DimNames = MRStruct.RecoPar.dimnames;
elseif(isfield(MRStruct.RecoPar,'dimnames_small'))
    DimNames = MRStruct.RecoPar.dimnames_small;
end




%%
InputName = inputname(1);
if(~isempty(InputName))
    evalin('caller',['clear ' InputName])
end



%% Perform Averaging

if(isfield(Settings,'AverageOverIndex'))
    AvgInd = Settings.AverageOverIndex;
elseif(exist('DimNames','var'))
    AvgInd =  find(~cellfun(@isempty,regexpi(DimNames,'Ave')));
    if(isempty(AvgInd))
        fprintf('\nWarning in op_AverageMRData: There seems to be no dimension that corresponds to averages (name ''Ave''). Don''t do anything.')       
        return;
    end
else
    fprintf('\nWarning in op_AverageMRData: I am not sure which dimension I should average. Don''t do anything.')
    return;
end
MRStruct.Data = sum(MRStruct.Data,AvgInd)/size(MRStruct.Data,AvgInd);
MRStruct.RecoPar.DataSize = size(MRStruct.Data);
MRStruct.RecoPar.nAverages = 1;
if(isfield(MRStruct.RecoPar,'dimnames'))
    MRStruct.RecoPar.dimnames(AvgInd) = [];
elseif(isfield(MRStruct.RecoPar,'dimnames_small'))
    MRStruct.RecoPar.dimnames_small(AvgInd) = [];
end


%%

if(Settings.UndoWeightedAveraging_flag && isfield(MRStruct,'WeightedAveragingMap'))
    Tmp = MRStruct.WeightedAveragingMap;
    Tmp(Tmp == 0) = 1;
    MRStruct.Data = MRStruct.Data ./ Tmp; 
end

%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);


