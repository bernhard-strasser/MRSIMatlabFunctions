function [MPRStruct] = io_ReadMPR(MPRStruct,Settings)
%
% op_ReadAndRecoBorjanSpiralData Read and reconstruct data from Borjan Gagoski's Spiral MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in Spiral MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         SpiralDataFile          ...     
% -         SpiralTrajectoryFile    ...     
% -         Settings                ...     
%
% Output:
% -         ?                      ...     
% -         ?                        ...     
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

% Input data must be structure with at least fields 'Data' and 'Par'.
% 'Data' must be of size 



%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end


%% Read MPR Data

files = dir(MPRStruct.DataFile);
FoldersToRead = files([files(:).isdir]);
FoldersToRead = FoldersToRead(~ismember({FoldersToRead(:).name},{'.','..'}));
FoldersToRead = {FoldersToRead.name};
if(isempty(FoldersToRead))
    FoldersToRead{1} = [];
end

% Read Parameters from any file
files = dir([MPRStruct.DataFile '/' FoldersToRead{1}]);
MPRStruct.Par = read_ascconv([MPRStruct.DataFile '/' FoldersToRead{1} '/' files(3).name]);
MPRStruct.Par.ContrastNames = FoldersToRead;
MPRStruct.Par.dimnames = {'x','y','z','t','cha','contrast'};
MPRStruct.RecoPar = MPRStruct.Par;      % The data was alrdy reconstructed by ICE


MPRStruct.Data = dicomreadVolume([MPRStruct.DataFile '/' FoldersToRead{1}]);
for CurFold = 2:numel(FoldersToRead)
     MPRStruct.Data(:,:,:,:,:,CurFold) = dicomreadVolume([MPRStruct.DataFile '/' FoldersToRead{CurFold}]);
end

MPRStruct.Data = double(permute(MPRStruct.Data, [2 4 1 3 5 6]));
MPRStruct.Data = MPRStruct.Data(:,:,end:-1:1,:,:,:,:);
% Data Size: x,y,z,t,cha,contrast



clear files FoldersToRead CurFold


%% Create Mask

INV2Logi = (strcmpi(MPRStruct.Par.ContrastNames,'INV2'));

if(any(INV2Logi))
    MPRStruct.Mask = MPRStruct.Data(:,:,:,:,:,INV2Logi) > max(MPRStruct.Data(:,:,:,:,:,INV2Logi))/25;
end
clear INV2Logi


%% Postparations


MPRStruct = supp_UpdateRecoSteps(MPRStruct,Settings);
MPRStruct = supp_FixPars(MPRStruct);



