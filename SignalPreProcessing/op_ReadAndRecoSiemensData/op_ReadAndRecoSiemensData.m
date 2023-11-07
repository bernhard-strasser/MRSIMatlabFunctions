function [MRStruct,RefStruct,NoiseStruct] = op_ReadAndRecoSiemensData(file,NonCartTrajFile,AdditionalInput,Settings)
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


%% 0. Preparations

if(~exist('Settings','var'))
    AdditionalInput = struct();
end
if(~exist('Settings','var'))
    Settings = struct();
end
if(~exist('NonCartTrajFile','var'))
    NonCartTrajFile = [];
end


% Always make cell to be consistent
if(~iscell(file))
    tmp{1} = file;
    file = tmp; clear tmp;
end

if(numel(file) > 1 && ~endsWith(file{1},'IMA'))
    fprintf('\nError in op_ReadAndRecoSiemensData: You gave me several non-DICOM files. I cannot digest that. Need to stop here.')
    return;
end


%% Read Data

dicom_flag = false;
if(  endsWith(file{1},'IMA') || isdir(file{1})  )
    dicom_flag = true;
    MRStruct = read_csi_dicom(file);
    MRStruct.Par.dicom_flag = true;	
    RefStruct = struct; NoiseStruct = struct;
else
    [MRStruct,RefStruct,NoiseStruct] = io_ReadAndReshapeSiemensData(file{1},NonCartTrajFile);
    MRStruct.Par.dicom_flag = false;
end



%% Handle Noisedata

if(~(isfield(NoiseStruct,'Data') && ~isempty(NoiseStruct.Data)) && isfield(AdditionalInput,'NoiseData'))
    NoiseStruct.Data = AdditionalInput.NoiseData; 
end


%% Reconstruct Data
if(~dicom_flag)
    MRStruct = op_ReconstructMRData(MRStruct,NoiseStruct,Settings);
    RefStruct = op_ReconstructMRData(RefStruct,NoiseStruct,Settings);
end

%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);
RefStruct = supp_UpdateRecoSteps(RefStruct,Settings);
NoiseStruct = supp_UpdateRecoSteps(NoiseStruct,Settings);



