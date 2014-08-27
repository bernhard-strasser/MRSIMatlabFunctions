function Noise_mat = GatherNoiseFromCSI(InData,Full_ElliptWeighted_Or_Weighted_Acq,nFIDendpoints,TakeNPointsOutOfEnd)
%
% read_csi_dat Read in csi-data from Siemens raw file format
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [csi,NoiseCorrMat,Noise_mat,InData] = read_csi_dat(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag, NoiseCorrMat)
%
% Input: 
% -         csi_path                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
% -         NoFFT_flag                  ...     If this is true, don't perform any fft.
% -         NoiseCorrMat                ...     If size(NoiseCorrMat) = [cha cha]: the k-space Data gets decorrelated with this matrix. 
%                                               If NoiseCorrMat = 1: the end of the FIDs of the outermost k-space/image points are used for noise decorrelation.
%                                               If NoiseCorrMat = 0, or not existant: No Noise Decorrelation performed
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         NoiseCorrMat                ...     The Noise Correlation Matrix in order to check if it was properly computed from the csi data. Is 0 if no decorrelation performed.
% -         Noise_mat                   ...     The Noise gathered from the CSI data. Noise_mat = 0, if not gathered.
% -         InData                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh, read_ascconv, hadamard_encoding.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


if(nargin < 1)
	fprintf('\nProblem: Cannot gather noise if no CSI or NoiseCorrMat is provided.')
	Noise_mat = 0;
	return;
end
if(~exist('Full_ElliptWeighted_Or_Weighted_Acq','var'))
	Full_ElliptWeighted_Or_Weighted_Acq = 2;	% = Elliptical
end
if(~exist('nFIDendpoints','var'))
	nFIDendpoints = 400;
end
if(~exist('TakeNPointsOutOfEnd','var'))
	TakeNPointsOutOfEnd = 200;	% = Elliptical
end


%% 1. Gather Noise


% Gather "Noise" of the CSI data (end of FIDs in the outermost acquired circle of k-space points). Only working for 2D and Multislice
randy = randperm(nFIDendpoints); randy = randy(1:TakeNPointsOutOfEnd); % Take 'em randomly

% Take random points at end of FID
Noise_csi = InData(:,:,:,:,:,end - (nFIDendpoints - 1) : end); Noise_csi = Noise_csi(:,:,:,:,:,randy);

% Take only csi-points which are farest away from k-space center (so a circle with radius 31 k-space points)
if(Full_ElliptWeighted_Or_Weighted_Acq == 2)
	[Elliptical_dummy,OuterkSpace_mask1] = EllipticalFilter(squeeze(InData(1,:,:,1,1,1,1)),[1 2],[1 1 1 size(InData,3)/2-1],1);
	[Elliptical_dummy,OuterkSpace_mask2] = EllipticalFilter(squeeze(InData(1,:,:,1,1,1,1)),[1 2],[1 1 1 size(InData,3)/2-2],1);
	OuterkSpace_mask = OuterkSpace_mask1 - OuterkSpace_mask2;
else
	OuterkSpace_mask = zeros([size(InData,2), size(InData,3)]);
	OuterkSpace_mask(1:end,1) = 1; OuterkSpace_mask(1,1:end) = 1; OuterkSpace_mask(end,1:end) = 1; OuterkSpace_mask(1:end,end) = 1;
end
PI_mask = abs(squeeze(InData(1,:,:,1,1,1))); PI_mask(PI_mask > 0) = 1;
csi_mask = OuterkSpace_mask .* PI_mask;
csi_mask = repmat(logical(csi_mask), [1 1 size(InData,1)*size(InData,4)*TakeNPointsOutOfEnd]);
csi_mask = reshape(csi_mask, [size(InData,2) size(InData,3) size(InData,1) size(InData,4) TakeNPointsOutOfEnd]);
csi_mask = permute(csi_mask, [3 1 2 4 5]);

Noise_mat = Noise_csi(csi_mask);
Noise_mat = reshape(Noise_mat, [size(InData,1) numel(Noise_mat)/size(InData,1)]);  




