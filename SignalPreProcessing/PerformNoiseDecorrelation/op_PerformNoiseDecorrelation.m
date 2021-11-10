function MRStruct = op_PerformNoiseDecorrelation(MRStruct,NoiseCorrMatStruct)
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
% [csi,NoiseCorrMat,Noise_mat,kSpace.CSI] = read_csi_dat(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag, NoiseCorrMat)
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
% -         kSpace.CSI                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


if(nargin < 2)
	fprintf('\nProblem: Cannot Noise-Decorrelate if no Input Data or NoiseCorrMat is provided.')
end

if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end


%% Find out where cha-dimension is
if(isfield(MRStruct.RecoPar,'dimnames_small'))
   ChaDim = find(strcmpi(MRStruct.RecoPar.dimnames_small,'cha'));
else    % Assume Cha is first dimension
    ChaDim = 1;                
end
if(ChaDim ~= 1 && ChaDim ~= numel(size(MRStruct.Data)))
    error('Error in op_PerformNoiseDecorrelation: Channel dimension is neither first nor last')
end

%% 1. Memory Considerations - Find best Looping

size_InData = size(MRStruct.Data);

[dummy, MemFree] = memused_linux(1);
MemNecessary = numel(MRStruct.Data)*8*3*2/2^20;						% every entry of InData is double --> 8 bytes. *2 because complex. *3 as safety measure (so InData fits 2 times in memory,
																% once it is already there and 2 more times it should fit in). /2^20 to get from bytes to MB.
if(MemNecessary > MemFree)
	LoopOverIndex = MemFree ./ (MemNecessary./size_InData);	% Because the 1st index is the coil. We can never loop over the coils.
	LoopOverIndex(LoopOverIndex < 1) = NaN; LoopOverIndex(ChaDim) = NaN;
	LoopOverIndex = find(nanmin(LoopOverIndex));
	LoopOverIndex = LoopOverIndex(1);								% Only take the 1st if several are the minimum.
	AccessString = [repmat(':,',[1 LoopOverIndex-1]) 'LoopIndex,' repmat(':,',[1 numel(size_InData)-LoopOverIndex])];
	AccessString(end) = [];
end


%% Scale NoiseCorrMat

NoiseScalingFactor = NoiseCorrMatStruct.Par.Dwelltimes(1) / MRStruct.RecoPar.Dwelltimes(1);
NoiseCorrMatStruct.Data = NoiseCorrMatStruct.Data * NoiseScalingFactor;


%% 1. Decorrelate

NoiseCorrMat_Chol = chol(NoiseCorrMatStruct.Data/2,'lower');

% Decorrelate within Loop due to memory usage
if(MemNecessary > MemFree)			
	
	%InData = zeros(size_InData);
	% Perform Noise Decorrelation Slice by Slice and Part by Part to avoid extensive memory usage 
	for LoopIndex = 1:size(MRStruct.Data,LoopOverIndex)

		TempData = eval(['MRStruct.Data(' AccessString ');']);
		size_TempData = size(TempData);
        
        
        if(ChaDim == 1)
            TempData = reshape(TempData, [size(TempData,1) numel(TempData)/size(TempData,1)]);
            TempData = NoiseCorrMat_Chol \ TempData;
        else
            TempData = reshape(TempData, [numel(TempData)/size(TempData,ChaDim) size(TempData,ChaDim)]);
            TempData = TempData / transpose(NoiseCorrMat_Chol);
        end
		eval(['MRStruct.Data(' AccessString ') = reshape(TempData,size_TempData);']);

	end
		
% Decorrelate everything at once.
else			
    if(ChaDim == 1)
        MRStruct.Data = reshape(MRStruct.Data, [size_InData(1) numel(MRStruct.Data)/size_InData(1)]);
        MRStruct.Data = NoiseCorrMat_Chol \ MRStruct.Data;
        MRStruct.Data = reshape(MRStruct.Data,size_InData);
    else
        MRStruct.Data = reshape(MRStruct.Data, [numel(MRStruct.Data)/size_InData(end) size_InData(end)]);
        MRStruct.Data = MRStruct.Data / transpose(NoiseCorrMat_Chol);
        MRStruct.Data = reshape(MRStruct.Data,size_InData);        
    end

end




