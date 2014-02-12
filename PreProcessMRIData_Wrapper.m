function [iSpace,Noise,PreProcessingInfo, kSpace] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo)
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
% [csi,NoiseCorrMat,Noise_mat,kSpace.CSI] = read_csi_dat(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, iSpaceShift,iSpaceShift(2),NoFFT_flag, NoiseCorrMat)
%
% Input: 
% -         csi_path                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         iSpaceShift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by iSpaceShift voxels
% -         iSpaceShift(2)                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by iSpaceShift(2) voxels
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
% File dependancy: memused_linux_1_0,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


% Find out memory used by MATLAB
memused_before = memused_linux_1_0(1); 



% % Assign standard values to variables if nothing is passed to function.
% if(~exist('zerofill_to_nextpow2_flag','var'))
%     zerofill_to_nextpow2_flag = 1;
% end
% if(~exist('zerofilling_fact','var'))
%     zerofilling_fact = 1;
% end
% if(~exist('Hadamard_flag','var'))
%     Hadamard_flag = 0;
% end
% if(~exist('iSpaceShift','var'))
%     iSpaceShift = 0;
% end
% if(~exist('iSpaceShift(2)','var'))
%     iSpaceShift(2) = 0;
% end
% if(~exist('NoFFT_flag','var'))
%     NoFFT_flag = false;
% end
% if(~exist('NoiseCorrMat','var'))
%     NoiseCorrMat = false;
% end
% if(~exist('Noise_mat','var'))
%     Noise_mat = 0;
% end



%% 1. Loop Over All Existing DataSets In PreProcessingInfo and Preprocess Data

DataSetNames = transpose(fields(PreProcessingInfo));
for CurDataSet = DataSetNames

	if(~isfield(kSpace,CurDataSet{:}))
		continue;
	end
	
	
	% 1.0 InLoopPreps
	% Make String out of Cell
	CurDataSetString = CurDataSet{:};
	% Set fredirshift
	if(strcmpi(CurDataSetString, 'PATREFANDIMASCAN'))
		PreProcessingInfo.PATREFANDIMASCAN.fredir_shift = round(-ReadInInfo.General.Ascconv.Pos_Sag(1)/ReadInInfo.General.Ascconv.nFreqEnc(1));
	end
	
	
	
	
	
	
	
	% 1.1 Compute Noise Correlation Matrix if necessary
	% Only read in data if NoiseCorrelationMatrix exists and is a 1
	if( isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') && numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) == 1 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 1 )


		%%%%%%% GET NOISE %%%%%%

		% Get data from Noise Prescan
		if(isfield(kSpace,'NOISEADJSCAN'))
			fprintf('\nNoise Decorrelating Using\t...\tNoise Prescan.')    
			Noise_mat = reshape(kSpace.NOISEADJSCAN,[size(kSpace.NOISEADJSCAN,1) numel(kSpace.NOISEADJSCAN)/size(kSpace.NOISEADJSCAN,1)]);
		% Or from end of FID of the outer kSpace
		elseif(isfield(kSpace,'ONLINE'))
			fprintf('\nNoise Decorrelating Using\t...\tONLINE Data Itself.')
			Noise_mat = GatherNoiseFromCSI(kSpace.ONLINE,ReadInInfo.General.Ascconv.Full_ElliptWeighted_Or_Weighted_Acq);
		end

		% Copy NoiseCorrMat
		for CurDataSet2 = DataSetNames
			if(strcmpi(CurDataSet2{:},'NOISEADJSCAN'))
				continue;
			end

			% Rescale Noise to be the same like ONLINE noise
			if( isfield(ReadInInfo,'NOISEADJSCAN') && isfield(ReadInInfo.NOISEADJSCAN,'Dwelltime') && isfield(ReadInInfo,CurDataSet2{:}) && isfield(ReadInInfo.(CurDataSet2{:}),'Dwelltime') )
				NoiseScalingFactor = sqrt(ReadInInfo.(CurDataSet2{:}).Dwelltime / ReadInInfo.NOISEADJSCAN.Dwelltime);
			else
				NoiseScalingFactor = 1;
			end
			Noise_mat = Noise_mat * NoiseScalingFactor;

			% Compute noise correlation matrix
			NoiseCorrMat = 1/(size(Noise_mat,2)) * (Noise_mat * Noise_mat');
			PreProcessingInfo.(CurDataSet2{:}).NoiseCorrMat = NoiseCorrMat;
			Noise.(CurDataSet2{:}) = Noise_mat;
		end

		% Delete Uneccessary stuff
		clear NoisePar Noise_mat NoiseCorrMat NoiseScalingFactor;
	end

	
	
	
	
	% 1.2. Perform Noise Decorrelation
	if(isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') && numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) > 1)
		kSpace.(CurDataSetString) = PerformNoiseDecorrelation(kSpace.(CurDataSetString),PreProcessingInfo.(CurDataSetString).NoiseCorrMat);
	end

	

	% 1.2 Hadamard Decoding
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hadamard_flag') && PreProcessingInfo.(CurDataSetString).Hadamard_flag && size(kSpace.(CurDataSetString),5) > 1)
		kSpace.(CurDataSetString) = hadamard_decoding(kSpace.(CurDataSetString),5);            
	end
	
	
	
	% 1.3. Flip in kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceAlong') && isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceWhileAccessing'))
		for FlipAlong = PreProcessingInfo.(CurDataSetString).FlipkSpaceAlong
			eval(['kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = flipdim(kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,FlipAlong);']);
			CircshiftVec = zeros([1 numel(size(kSpace.(CurDataSetString)))]);
			CircshiftVec(FlipAlong) = 1; %#ok
			eval(['kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = circshift(kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,CircshiftVec);']);
		end
	end



	% 1.3. Circshift in kSpace
	if(	isfield(PreProcessingInfo.(CurDataSetString),'ShiftkSpace') && sum(~(PreProcessingInfo.(CurDataSetString).ShiftkSpace == 0)) > 0 && numel(PreProcessingInfo.(CurDataSetString).ShiftkSpace) == numel(size(kSpace.(CurDataSetString))))
		kSpace.(CurDataSetString) = circshift(kSpace.(CurDataSetString),PreProcessingInfo.(CurDataSetString).ShiftkSpace);		
	end





	
	% 1.5. Correct phase due to shift in freq encoding direction
	% Strange phenomenon for images: if you shift them in frequency encoding direction, the phase of the image changes. This undoes the phase change.
	if(isfield(PreProcessingInfo.(CurDataSetString),'fredir_shift') && ne(PreProcessingInfo.(CurDataSetString).fredir_shift(1),0))                                                    
		kSpace.(CurDataSetString) = kSpace.(CurDataSetString) * ...
		exp(1i*pi/180*PreProcessingInfo.(CurDataSetString).fredir_shift(1)*(360/2^nextpow2(size(kSpace.(CurDataSetString),2))));     % Is it really 2^nextpow2(fredir_measured) 
	end																																% and not just fredir_measured?                                                                                          

	
	
	% 1.5. Save Unfiltered kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'SaveUnfilteredkSpace') && PreProcessingInfo.(CurDataSetString).SaveUnfilteredkSpace)
		kSpace.([CurDataSetString '_Unfiltered']) = kSpace.(CurDataSetString);
		PreProcessingInfo.([CurDataSetString '_Unfiltered']) = PreProcessingInfo.(CurDataSetString);
	end

	
	% 8. Apply Elliptical Filter
	if(isfield(PreProcessingInfo.(CurDataSetString),'EllipticalFilterSize') && PreProcessingInfo.(CurDataSetString).EllipticalFilterSize > 0)
		OversamplingFactor = 1;				% Hack for now
		kSpace.(CurDataSetString) = EllipticalFilter_1_1(kSpace.(CurDataSetString),[2 3],[OversamplingFactor 1 1 PreProcessingInfo.(CurDataSetString).EllipticalFilterSize],1);
	end




	% 9. Apply Hamming Filter
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hamming_flag') && PreProcessingInfo.(CurDataSetString).Hamming_flag)
		kSpace.(CurDataSetString) = HammingFilter_1_3(kSpace.(CurDataSetString),[2 3],1);
	end

	
	
	
	

	% 1.3 FFT FROM K_SPACE TO DIRECT SPACE
	if(isfield(PreProcessingInfo.(CurDataSetString), 'NoFFT_flag') && ~PreProcessingInfo.(CurDataSetString).NoFFT_flag)
		
		if(strcmpi(CurDataSetString,'ONLINE'))
			ConjFlag = true;
		else
			ConjFlag = false;
		end
		iSpace.(CurDataSetString) = FFTOfMRIData(kSpace.(CurDataSetString),ConjFlag);
	else
		iSpace.(CurDataSetString) = 0;
	end




	% 6. Spatial Shift
	if(isfield(PreProcessingInfo.(CurDataSetString), 'iSpaceShift') && ne(PreProcessingInfo.(CurDataSetString).iSpaceShift(1),0))
		iSpace.(CurDataSetString) = circshift(iSpace.(CurDataSetString), [0 iSpaceShift(1) 0 0 0 0 0]);
	end
	if(isfield(PreProcessingInfo.(CurDataSetString), 'iSpaceShift') && numel(PreProcessingInfo.(CurDataSetString).iSpaceShift) > 1 && ne(PreProcessingInfo.(CurDataSetString).iSpaceShift(2),0))
		iSpace.(CurDataSetString) = circshift(iSpace.(CurDataSetString), [0 0 iSpaceShift(2) 0 0 0 0]);
	end



end



%% 7. Postparations

if(nargout > 1 && ~exist('Noise','var'))
	Noise = 0;
end
if(nargout > 2 && ~exist('kSpace','var'))
	kSpace = 0;
end
if(nargout > 3 && ~exist('PreProcessingInfo','var'))
	PreProcessingInfo = 0;
end

memused_after = memused_linux_1_0(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


