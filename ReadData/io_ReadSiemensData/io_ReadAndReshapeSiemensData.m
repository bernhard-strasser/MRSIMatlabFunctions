function [MRStruct,RefStruct,NoiseStruct] = io_ReadAndReshapeSiemensData(file,NonCartTrajFile)
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

RefStruct = struct();
NoiseStruct = struct();
Settings = struct();

if(~exist('NonCartTrajFile','var'))
    NonCartTrajFile = [];
end


% Always make cell to be consistent
if(~iscell(file))
    tmp{1} = file;
    file = tmp; clear tmp;
end


%% Read Data

if(endsWith(file{1},'IMA') || isdir(file{1}))
	MRStruct = read_csi_dicom(file); % Differences between VB and VD/VE?
	MRStruct.Par.dicom_flag = true;	
    RefStruct = struct; NoiseStruct = struct;
else
    MRStruct.Par = read_ascconv(file{1});
    MRStruct.Par.dicom_flag = false;	
	if(strcmpi(MRStruct.Par.AssumedSequence,'ViennaCRT'))
		[MRStruct,RefStruct,NoiseStruct] = io_ReadAverageReshape3DCRTDataOwnRead(file{1});
        MRStruct.TrajFile = NonCartTrajFile; RefStruct.TrajFile = NonCartTrajFile;
	elseif(strcmpi(MRStruct.Par.AssumedSequence,'BorjanSpiral'))
		MRStruct = io_ReadAverageReshapeBorjanSpiralData(file{1},NonCartTrajFile);
    elseif(strcmpi(MRStruct.Par.AssumedSequence,'CSIOrSVS'))
        [MRStruct,RefStruct,NoiseStruct] = io_ReadAverageReshapePhaseEncodedMRSIData(file{1});
    elseif(strcmpi(MRStruct.Par.AssumedSequence,'Imaging_GRE'))
        MRStruct = io_ReadAverageReshapeImagingGREData(file{1});
    elseif(strcmpi(MRStruct.Par.AssumedSequence,'AntoinesEccentricOrRosette')) % SB2022
        MRStruct = io_ReadAverageReshape3DEccentricDataOwnRead(file{1});   
	else
		MRStruct = io_ReadAverageReshapeGenericData(file{1});
	end
end


%% Add files

MRStruct.DataFile = file;
if(exist('NonCartTrajFile','var'))
    MRStruct.TrajFile = NonCartTrajFile;
end
if(~isempty(RefStruct) && ~isempty(fieldnames(RefStruct)))
    RefStruct.DataFile = file;
    if(exist('NonCartTrajFile','var'))
        RefStruct.TrajFile = NonCartTrajFile;
    end
end
if(~isempty(NoiseStruct) && ~isempty(fieldnames(NoiseStruct)))
    NoiseStruct.DataFile = file;
end

%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);
RefStruct = supp_UpdateRecoSteps(RefStruct,Settings);
NoiseStruct = supp_UpdateRecoSteps(NoiseStruct,Settings);



