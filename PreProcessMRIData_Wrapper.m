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




%% 1. Loop Over All Existing DataSets In PreProcessingInfo and Preprocess Data

DataSetNames = transpose(fields(PreProcessingInfo));
for CurDataSet = DataSetNames

	if(~isfield(kSpace,CurDataSet{:}))
		continue;
	end
	fprintf('\n\nPreprocessing Data         \t...\t%s',CurDataSet{:})
	
	
	% 1.0 InLoopPreps
	% Make String out of Cell
	CurDataSetString = CurDataSet{:};
	% Set fredirshift
	if(~isfield(PreProcessingInfo.(CurDataSetString),'fredir_shift') && size(kSpace.(CurDataSetString),6) == 1)	% Only for imaging data
		PreProcessingInfo.(CurDataSetString).fredir_shift = 2*ReadInInfo.General.Ascconv.Pos_Sag(1)/ (-ReadInInfo.General.Ascconv.FoV_Read(1)/ReadInInfo.(CurDataSetString).nReadEnc );
	end
	
	
	% Do conj if necessary
	
	
	
	
	% 1.1 Compute Noise Correlation Matrix if necessary
	% Only read in data if NoiseCorrelationMatrix exists and is a 1
	if( isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') && numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) == 1 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat > 0 )
		fprintf('\nCompute noise correlation matrix')


		%%%%%%% GET NOISE %%%%%%

		% Get data from Noise Prescan
		if(isfield(kSpace,'NOISEADJSCAN') && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 1)
			fprintf('\nGet noise from\t...\tNoise Prescan.')    
			Noise_mat = reshape(kSpace.NOISEADJSCAN,[size(kSpace.NOISEADJSCAN,1) numel(kSpace.NOISEADJSCAN)/size(kSpace.NOISEADJSCAN,1)]);
			Noise_mat = Noise_mat(:,20:end);
		% Or from end of FID of the outer kSpace
		elseif(isfield(kSpace,'ONLINE') && size(kSpace.ONLINE,6) > 512 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 2)
			fprintf('\nGet noise from\t...\tONLINE Data Itself.')
			Noise_mat = GatherNoiseFromCSI(kSpace.ONLINE,ReadInInfo.General.Ascconv.Full_ElliptWeighted_Or_Weighted_Acq);
		end

		% Copy NoiseCorrMat
		if(exist('Noise_mat','var'))
			for CurDataSet2 = DataSetNames
				if(strcmpi(CurDataSet2{:},'NOISEADJSCAN'))
					continue;
				end

				% Rescale Noise to be the same like ONLINE noise
				if( numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) == 1 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 1 && ...
					isfield(ReadInInfo,'NOISEADJSCAN') && isfield(ReadInInfo.NOISEADJSCAN,'Dwelltime') && isfield(ReadInInfo,CurDataSet2{:}) && isfield(ReadInInfo.(CurDataSet2{:}),'Dwelltime') )
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
	end


	
	% 1.2. Perform Noise Decorrelation
	if(isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') && numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) > 1)
		fprintf('\nPerform Noise Decorrelation.')    		
		kSpace.(CurDataSetString) = PerformNoiseDecorrelation(kSpace.(CurDataSetString),PreProcessingInfo.(CurDataSetString).NoiseCorrMat);
	end

	

	% 1.2 Hadamard Decoding
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hadamard_flag') && PreProcessingInfo.(CurDataSetString).Hadamard_flag && size(kSpace.(CurDataSetString),5) > 1)
		kSpace.(CurDataSetString) = hadamard_decoding(kSpace.(CurDataSetString),5);            
	end
	
	
	
	% 1.3. Flip in kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceAlong') && isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceWhileAccessing'))
		fprintf('\nFlip in kSpace.')    		
		for FlipAlong = PreProcessingInfo.(CurDataSetString).FlipkSpaceAlong
			try
				eval(['kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = flipdim(kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,FlipAlong);']);
			catch errie
				fprintf('\nProblem in PreProcessMRIData_Wrapper: Tried to execute\n[kSpace.%s(%s) = flipdim(kSpace.%s(%s) ,FlipAlong);\nThis resulted in an error:\n''%s'' at line %d.\n', ...
				CurDataSetString, PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing,CurDataSetString,...
				PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing,errie.message,errie.stack(1).line);
			end
			if(~exist('errie','var'))
				CircshiftVec = zeros([1 numel(size(kSpace.(CurDataSetString)))]);
				CircshiftVec(FlipAlong) = 1; %#ok
				eval(['kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = circshift(kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,CircshiftVec);']);
			end
			clear errie;
		end
	end



	% 1.3. Circshift in kSpace
	if(	isfield(PreProcessingInfo.(CurDataSetString),'ShiftkSpace') && sum(~(PreProcessingInfo.(CurDataSetString).ShiftkSpace == 0)) > 0 && numel(PreProcessingInfo.(CurDataSetString).ShiftkSpace) == numel(size(kSpace.(CurDataSetString))))
		fprintf('\nCircshift in kSpace.')    				
		kSpace.(CurDataSetString) = circshift(kSpace.(CurDataSetString),PreProcessingInfo.(CurDataSetString).ShiftkSpace);		
	end


	% 1.4 Correct global phase
	if(isfield(PreProcessingInfo.(CurDataSetString),'GlobalPhase_rad') && ne(PreProcessingInfo.(CurDataSetString).GlobalPhase_rad(1),0)) 
		fprintf('\nAdd global phase of %f in kSpace.',rad2deg(PreProcessingInfo.(CurDataSetString).GlobalPhase_rad(1)))    						
		kSpace.(CurDataSetString) = kSpace.(CurDataSetString) * exp(1i*PreProcessingInfo.(CurDataSetString).GlobalPhase_rad(1));
	end		


	% 1.5. Correct phase due to shift in freq encoding direction
	% Strange phenomenon for images: if you shift them in frequency encoding direction, the phase of the image changes. This undoes the phase change.
	if(isfield(PreProcessingInfo.(CurDataSetString),'fredir_shift') && ne(PreProcessingInfo.(CurDataSetString).fredir_shift(1),0))
		phase = PreProcessingInfo.(CurDataSetString).fredir_shift(1)*(360/2^nextpow2(size(kSpace.(CurDataSetString),2)));
		fprintf('\nCorrect global phase of %f deg due to FoV-shift of %f in kSpace.',phase,PreProcessingInfo.(CurDataSetString).fredir_shift(1)) 								
		kSpace.(CurDataSetString) = kSpace.(CurDataSetString) * exp(1i*deg2rad(phase));												% Is it really 2^nextpow2(fredir_measured) 
		clear phase;																												% and not just fredir_measured?  
	end                                                                                        

	
	
	% 1.5. Save Unfiltered kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'SaveUnfilteredkSpace') && PreProcessingInfo.(CurDataSetString).SaveUnfilteredkSpace)
		fprintf('\nSave Unfiltered kSpace.')    				
		kSpace.([CurDataSetString '_Unfiltered']) = kSpace.(CurDataSetString);
		PreProcessingInfo.([CurDataSetString '_Unfiltered']) = PreProcessingInfo.(CurDataSetString);
	end

	
	
	% 8. Apply Elliptical Filter
	if(isfield(PreProcessingInfo.(CurDataSetString),'EllipticalFilterSize') && PreProcessingInfo.(CurDataSetString).EllipticalFilterSize(end) > 0)
		fprintf('\nApply Elliptical Filter with Filtersize %d.',transpose(PreProcessingInfo.(CurDataSetString).EllipticalFilterSize))    				
		kSpace.(CurDataSetString) = EllipticalFilter_1_1(kSpace.(CurDataSetString),[2 3],PreProcessingInfo.(CurDataSetString).EllipticalFilterSize,1);
	end




	% 9. Apply Hamming Filter
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hamming_flag') && PreProcessingInfo.(CurDataSetString).Hamming_flag)
		fprintf('\nApply Hamming Filter.')    				
		kSpace.(CurDataSetString) = HammingFilter(kSpace.(CurDataSetString),[2 3],1);
	end

	
	% Zerofilling in kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'ZeroFillingDesiredSize') && numel(PreProcessingInfo.(CurDataSetString).ZeroFillingDesiredSize) > 1)
		fprintf('\nZerofill in kSpace from size [%s] to [%s]',sprintf('%d ',size(kSpace.(CurDataSetString))),sprintf('%d ',PreProcessingInfo.(CurDataSetString).ZeroFillingDesiredSize))    						
		kSpace.(CurDataSetString) = Zerofilling(kSpace.(CurDataSetString), PreProcessingInfo.(CurDataSetString).ZeroFillingDesiredSize);
	end
	

	% 1.3 FFT FROM K_SPACE TO DIRECT SPACE
	if(isfield(PreProcessingInfo.(CurDataSetString), 'NoFFT_flag') && ~PreProcessingInfo.(CurDataSetString).NoFFT_flag)
		fprintf('\nFFT Data.')    						
		if(strcmpi(CurDataSetString,'ONLINE') && size(kSpace.(CurDataSetString),6) > 1)
			ConjFlag = true;
		else
			ConjFlag = false;
		end
		iSpace.(CurDataSetString) = FFTOfMRIData(kSpace.(CurDataSetString),ConjFlag);
	else
		iSpace.(CurDataSetString) = 0;
	end


	
	% Remove spatial oversampling
	if(isfield(PreProcessingInfo.(CurDataSetString), 'RmOs') && PreProcessingInfo.(CurDataSetString).RmOs == 1)
		fprintf('\nRemove Spatial Oversampling.')
		image_center = floor(size(iSpace.(CurDataSetString),2) / 2) + 1;
		left_border = image_center - ceil(size(iSpace.(CurDataSetString),2) / (2*2));   
		right_border = image_center + ceil(size(iSpace.(CurDataSetString),2) / (2*2)) - 1;   
	
		if(size(iSpace.(CurDataSetString),2) > 1)
			iSpace.(CurDataSetString) = iSpace.(CurDataSetString)(:,left_border:right_border,:,:,:,:,:);
		else
			iSpace.(CurDataSetString) = fftshift(fft(ifftshift(kSpace.(CurDataSetString),2),[],2),2);
			iSpace.(CurDataSetString) = iSpace.(CurDataSetString)(:,left_border:right_border,:,:,:,:,:);
		end
		if(nargout > 3)
			kSpace.(CurDataSetString) = ifftshift(ifft(fftshift(iSpace.(CurDataSetString),2),[],2),2);
			if(isfield(kSpace,[CurDataSetString '_Unfiltered']))
				kSpace.([CurDataSetString '_Unfiltered']) = kSpace.([CurDataSetString '_Unfiltered'])(:,1:2:end,:,:,:,:,:);
			end
		end	
	end
	


	% 6. Spatial Shift
	if(isfield(PreProcessingInfo.(CurDataSetString), 'iSpaceShift') && ne(PreProcessingInfo.(CurDataSetString).iSpaceShift(1),0))
		fprintf('\nSpatial shift data by [%s].', sprintf('%d ',[0 iSpaceShift(1) 0 0 0 0 0]))    								
		iSpace.(CurDataSetString) = circshift(iSpace.(CurDataSetString), [0 iSpaceShift(1) 0 0 0 0 0]);
	end
	if(isfield(PreProcessingInfo.(CurDataSetString), 'iSpaceShift') && numel(PreProcessingInfo.(CurDataSetString).iSpaceShift) > 1 && ne(PreProcessingInfo.(CurDataSetString).iSpaceShift(2),0))
		fprintf('\nSpatial shift data by [%s].', sprintf('%d ',[0 0 iSpaceShift(2) 0 0 0 0]))    								
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


