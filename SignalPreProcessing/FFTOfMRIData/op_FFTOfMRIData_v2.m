function MRStruct = op_FFTOfMRIData_v2(MRStruct,Settings)
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
% InData = FFTOfMRIData(InData,Settings.ConjFlag,Settings.ApplyAlongDims,Settings.Ifft_flag,Settings.quiet_flag,Settings.FlipDim_flag)
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
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations

if(nargin < 1)
	fprintf('\nProblem: Cannot fft non-existant data. Provide some data.')
	MRStruct = struct;
	return;
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

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'ConjFlag'))
	Settings.ConjFlag = false;
end
if(~isfield(Settings,'ApplyAlongDims'))
	Settings.ApplyAlongDims = 2:min(numel(MRStruct.Par.DataSize),4);    % If enough dimensions are available, use 2:4
end
if(~isfield(Settings,'Ifft_flag'))
	Settings.Ifft_flag = false;
end
if(~isfield(Settings,'quiet_flag'))
    Settings.quiet_flag = false;
end
if(~isfield(Settings,'FlipDim'))
    Settings.FlipDim = 1;
end


%% Early Exit 
% If size is degenerate along dims along which FFT should be applied 

if(all(MRStruct.Par.DataSize(Settings.ApplyAlongDims) == 1))
    if(Settings.ConjFlag)
        MRStruct.Data = conj(MRStruct.Data);
        if(isfield(MRStruct,'NoiseData'))
            MRStruct.NoiseData = conj(MRStruct.NoiseData);
        end
        fprintf('\nDid only conj, no fft.')
    end
	return;
end


%% Call FFT Function

MRStruct.Data = FFTOfMRIData(MRStruct.Data,Settings.ConjFlag,Settings.ApplyAlongDims,Settings.Ifft_flag,Settings.quiet_flag,Settings.FlipDim);
if(isfield(MRStruct,'NoiseData'))
    MRStruct.NoiseData = FFTOfMRIData(MRStruct.NoiseData,Settings.ConjFlag,Settings.ApplyAlongDims,Settings.Ifft_flag,Settings.quiet_flag,Settings.FlipDim);    
end


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



