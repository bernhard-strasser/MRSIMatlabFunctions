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
% [iSpace,Noise,PreProcessingInfo, kSpace] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo)
%
% Input: 
% -         kSpace                      ...     Path of MRS(I) file.
% -         PreProcessingInfo           ...     Info about how the read in data should be pre-processed.
%                                               Sub-Fields: For each subdataset (EvalInfoMask entry) one, e.g. 'ONLINE', 'PATREFANDIMASCAN', and 'NOISEADJSCAN'
%                                               Sub-Sub-Fields: Each Sub-Field can have the following entries:
%                                               - NoFFT_flag: If true, do no fft is performed from kSpace to iSpace. Example: PreProcessingInfo.ONLINE.NoFFT_flag = true;
%                                               - fredir_shift: Corrects for a different phase caused by a shift in the frequency encoding direction
%                                               - FlipkSpaceAlong = 2;
%												- FlipkSpaceWhileAccessing = ':,:,:,:,:,:,2';
%                                               - (tbc)
% -         ReadInInfo                  ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% Output:
% -         iSpace                      ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         Noise                       ...     The Noise Correlation Matrix in order to check if it was properly computed from the csi data. Is 0 if no decorrelation performed.
% -         PreProcessingInfo           ...     The  updated PreProcessingInfo, as it is changed when reading in.
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh, read_ascconv, hadamard_encoding.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations


% Find out memory used by MATLAB
memused_before = memused_linux(1); 


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
	if(~isfield(PreProcessingInfo.(CurDataSetString),'fredir_shift') && size(kSpace.(CurDataSetString){1},6) == 1)	% Only for imaging data
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
			Noise_mat = reshape(kSpace.NOISEADJSCAN{1},[size(kSpace.NOISEADJSCAN{1},1) numel(kSpace.NOISEADJSCAN{1})/size(kSpace.NOISEADJSCAN{1},1)]);
			Noise_mat = Noise_mat(:,20:end);
		% Or from end of FID of the outer kSpace
		elseif(isfield(kSpace,'ONLINE') && size(kSpace.ONLINE{1},6) > 512 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 2)
			fprintf('\nGet noise from\t...\tONLINE Data Itself.')
			Noise_mat = GatherNoiseFromCSI(kSpace.ONLINE{1},ReadInInfo.General.Ascconv.Full_ElliptWeighted_Or_Weighted_Acq);
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
					NoiseScalingFactor = sqrt(ReadInInfo.NOISEADJSCAN.Dwelltime / ReadInInfo.(CurDataSet2{:}).Dwelltime);
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
		for echo = 1:numel(kSpace.(CurDataSetString))
			kSpace.(CurDataSetString){echo} = PerformNoiseDecorrelation(kSpace.(CurDataSetString){echo},PreProcessingInfo.(CurDataSetString).NoiseCorrMat);
		end
	end

	

	% 1.2 Hadamard Decoding
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hadamard_flag') && PreProcessingInfo.(CurDataSetString).Hadamard_flag && size(kSpace.(CurDataSetString){1},5) > 1)
		for echo = 1:numel(kSpace.(CurDataSetString))		
			kSpace.(CurDataSetString){echo} = hadamard_decoding(kSpace.(CurDataSetString){echo},5);    
		end
	end
	
	
	
	% 1.3. Circshift in kSpace
	if(	isfield(PreProcessingInfo.(CurDataSetString),'ShiftkSpace') && sum(~(PreProcessingInfo.(CurDataSetString).ShiftkSpace == 0)) > 0)
		if(numel(PreProcessingInfo.(CurDataSetString).ShiftkSpace) == numel(size(kSpace.(CurDataSetString){1})))
			fprintf('\nCircshift in kSpace.')   
			for echo = 1:numel(kSpace.(CurDataSetString))		
				kSpace.(CurDataSetString){echo} = circshift(kSpace.(CurDataSetString){echo},PreProcessingInfo.(CurDataSetString).ShiftkSpace);		
			end
		else
			fprintf('\nError: Cannot shift in kSpace, please specify PreProcessingInfo.%s.ShiftkSpace with numel(...) = %d',CurDataSetString,numel(size(kSpace.(CurDataSetString){1})))
		end
	end
	
	
	
	% 1.4. Flip in kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceAlong') && isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceWhileAccessing'))
		fprintf('\nFlip in kSpace.')
		for echo = 1:numel(kSpace.(CurDataSetString))
			for FlipAlong = PreProcessingInfo.(CurDataSetString).FlipkSpaceAlong
				try
					eval(['kSpace.(CurDataSetString){echo}(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = flipdim(kSpace.(CurDataSetString){echo}(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,FlipAlong);']);
				catch errie
					fprintf('\nProblem in PreProcessMRIData_Wrapper: Tried to execute\n[kSpace.%s(%s) = flipdim(kSpace.%s(%s) ,FlipAlong);\nThis resulted in an error:\n''%s'' at line %d.\n', ...
					CurDataSetString, PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing,CurDataSetString,...
					PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing,errie.message,errie.stack(1).line);
				end
				if(~exist('errie','var'))
					CircshiftVec = zeros([1 numel(size(kSpace.(CurDataSetString){echo}))]);
					CircshiftVec(FlipAlong) = 1; %#ok
					eval(['kSpace.(CurDataSetString){echo}(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = circshift(kSpace.(CurDataSetString){echo}(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,CircshiftVec);']);
				end
				clear errie;
			end
		end
	end


	% 1.3. Circshift in kSpace
	if(	isfield(PreProcessingInfo.(CurDataSetString),'ShiftkSpacePostFlip') && sum(~(PreProcessingInfo.(CurDataSetString).ShiftkSpace == 0)) > 0)
		if(numel(PreProcessingInfo.(CurDataSetString).ShiftkSpace) == numel(size(kSpace.(CurDataSetString){1})))
			fprintf('\nCircshift in kSpace.')   
			for echo = 1:numel(kSpace.(CurDataSetString))		
				kSpace.(CurDataSetString){echo} = circshift(kSpace.(CurDataSetString){echo},PreProcessingInfo.(CurDataSetString).ShiftkSpace);		
			end
		else
			fprintf('\nError: Cannot shift in kSpace, please specify PreProcessingInfo.%s.ShiftkSpace with numel(...) = %d',CurDataSetString,numel(size(kSpace.(CurDataSetString){1})))
		end
	end



	% 1.4 Correct global phase
	if(isfield(PreProcessingInfo.(CurDataSetString),'GlobalPhase_rad') && ne(PreProcessingInfo.(CurDataSetString).GlobalPhase_rad(1),0)) 
		fprintf('\nAdd global phase of %f in kSpace.',rad2deg(PreProcessingInfo.(CurDataSetString).GlobalPhase_rad(1)))    
		for echo = 1:numel(kSpace.(CurDataSetString))
			kSpace.(CurDataSetString){echo} = kSpace.(CurDataSetString){echo} * exp(1i*PreProcessingInfo.(CurDataSetString).GlobalPhase_rad(1));
		end
	end		


	% 1.5. Correct phase due to shift in freq encoding direction
	% Strange phenomenon for images: if you shift them in frequency encoding direction, the phase of the image changes. This undoes the phase change.
	if(isfield(PreProcessingInfo.(CurDataSetString),'fredir_shift') && ne(PreProcessingInfo.(CurDataSetString).fredir_shift(1),0))
		phase = PreProcessingInfo.(CurDataSetString).fredir_shift(1)*(360/2^nextpow2(size(kSpace.(CurDataSetString),2)));
		fprintf('\nCorrect global phase of %f deg due to FoV-shift of %f in kSpace.',phase,PreProcessingInfo.(CurDataSetString).fredir_shift(1)) 	
		for echo = 1:numel(kSpace.(CurDataSetString))
			kSpace.(CurDataSetString){echo} = kSpace.(CurDataSetString){echo} * exp(1i*deg2rad(phase));												% Is it really 2^nextpow2(fredir_measured) 
			clear phase;																												% and not just fredir_measured?  
		end
	end                                                                                        

	
	
	% 1.5. Save Unfiltered kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'SaveUnfilteredkSpace') && PreProcessingInfo.(CurDataSetString).SaveUnfilteredkSpace)
		fprintf('\nSave Unfiltered kSpace.')    				
		kSpace.([CurDataSetString '_Unfiltered']) = kSpace.(CurDataSetString);
		PreProcessingInfo.([CurDataSetString '_Unfiltered']) = PreProcessingInfo.(CurDataSetString);
	end

	
	if(isfield(ReadInInfo,'General') && isfield(ReadInInfo.General,'Ascconv') && isfield (ReadInInfo.General.Ascconv,'AsymmetricEcho') && ReadInInfo.General.Ascconv.AsymmetricEcho && strcmpi(CurDataSetString,'ONLINE'))
		for echo = 1:numel(kSpace.(CurDataSetString))
			sizzy = size(kSpace.(CurDataSetString){echo});
			kSpace.(CurDataSetString){echo} = cat(2,zeros([sizzy(1) 2^nextpow2(sizzy(2))-sizzy(2) sizzy(3:end)]),kSpace.(CurDataSetString){echo});
		end
	end
	
	
	% 8. Apply Elliptical Filter
	if(isfield(PreProcessingInfo.(CurDataSetString),'EllipticalFilterSize') && PreProcessingInfo.(CurDataSetString).EllipticalFilterSize(end) > 0)
		fprintf('\nApply Elliptical Filter with Filtercoefficients [a b c R] = [%s].',num2str(PreProcessingInfo.(CurDataSetString).EllipticalFilterSize))
		for echo = 1:numel(kSpace.(CurDataSetString))
			kSpace.(CurDataSetString){echo} = EllipticalFilter(kSpace.(CurDataSetString){echo},[2 3],PreProcessingInfo.(CurDataSetString).EllipticalFilterSize,1);
		end
	end




	% 9. Apply Hamming Filter
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hamming_flag') && PreProcessingInfo.(CurDataSetString).Hamming_flag)
		fprintf('\nApply Hamming Filter.')
		if(isfield(PreProcessingInfo.(CurDataSetString),'HammingFactor'))
			FilterWidth = PreProcessingInfo.(CurDataSetString).HammingFactor;
		else
			FilterWidth = 100;
		end
		for echo = 1:numel(kSpace.(CurDataSetString))
			kSpace.(CurDataSetString){echo} = HammingFilter(kSpace.(CurDataSetString){echo},[2 3],FilterWidth,'Radial',1);
		end
	end
	
	% Zerofilling in kSpace
	if(isfield(PreProcessingInfo.(CurDataSetString),'ZeroFillingDesiredSize') && numel(PreProcessingInfo.(CurDataSetString).ZeroFillingDesiredSize{1}) > 1)
		for echo = 1:numel(kSpace.(CurDataSetString))
			fprintf('\nEcho %d: Zerofill in kSpace from size [%s] to [%s]',echo,sprintf('%d ',size(kSpace.(CurDataSetString){echo})),sprintf('%d ',PreProcessingInfo.(CurDataSetString).ZeroFillingDesiredSize{echo}))    						
			kSpace.(CurDataSetString){echo} = ZerofillOrCutkSpace(kSpace.(CurDataSetString){echo}, PreProcessingInfo.(CurDataSetString).ZeroFillingDesiredSize{echo});
		end
	end
	

	% 1.3 FFT FROM K_SPACE TO DIRECT SPACE
	if(  (isfield(PreProcessingInfo.(CurDataSetString), 'NoFFT_flag') && ~PreProcessingInfo.(CurDataSetString).NoFFT_flag) || ...
		(isfield(PreProcessingInfo.(CurDataSetString), 'RmOs') && PreProcessingInfo.(CurDataSetString).RmOs == 1)  )
		fprintf('\nFFT Data.')    						
		if(strcmpi(CurDataSetString,'ONLINE') && size(kSpace.(CurDataSetString){1},6) > 1)
			ConjFlag = true;
		else
			ConjFlag = false;
		end
		iSpace.(CurDataSetString) = cell([1 numel(kSpace.(CurDataSetString))]);
		for echo = 1:numel(kSpace.(CurDataSetString))
			iSpace.(CurDataSetString){echo} = FFTOfMRIData(kSpace.(CurDataSetString){echo},ConjFlag);
		end
    else
        if(nargout >= 4)
            iSpace.(CurDataSetString){1} = 0;            
        else
            iSpace = kSpace;
        end
	end


	
	% Remove spatial oversampling
	if(isfield(PreProcessingInfo.(CurDataSetString), 'RmOs') && PreProcessingInfo.(CurDataSetString).RmOs == 1)
		fprintf('\nRemove Spatial Oversampling.')
		for echo = 1:numel(kSpace.(CurDataSetString))
			image_center = floor(size(iSpace.(CurDataSetString){echo},2) / 2) + 1;
			left_border = image_center - ceil(size(iSpace.(CurDataSetString){echo},2) / (2*2));   
			right_border = image_center + ceil(size(iSpace.(CurDataSetString){echo},2) / (2*2)) - 1;   
			if(size(iSpace.(CurDataSetString){echo},2) > 1)
				iSpace.(CurDataSetString){echo} = iSpace.(CurDataSetString){echo}(:,left_border:right_border,:,:,:,:,:);
			else
				iSpace.(CurDataSetString){echo} = fftshift(fft(ifftshift(kSpace.(CurDataSetString){echo},2),[],2),2);
				iSpace.(CurDataSetString){echo} = iSpace.(CurDataSetString){echo}(:,left_border:right_border,:,:,:,:,:);
			end
			if(nargout > 3)
				kSpace.(CurDataSetString){echo} = FFTOfMRIData(iSpace.(CurDataSetString){echo},ConjFlag, [2 3 4],1);
				if(isfield(kSpace,[CurDataSetString '_Unfiltered']))
					kSpace.([CurDataSetString '_Unfiltered']){echo} = kSpace.([CurDataSetString '_Unfiltered']){echo}(:,1:2:end,:,:,:,:,:);
				end
			end
		end
	end
	


	% 6. Spatial Shift
	if(isfield(PreProcessingInfo.(CurDataSetString), 'iSpaceShift') && ne(PreProcessingInfo.(CurDataSetString).iSpaceShift(1),0))
		fprintf('\nSpatial shift data by [%s].', sprintf('%d ',[0 iSpaceShift(1) 0 0 0 0 0]))    		
		for echo = 1:numel(kSpace.(CurDataSetString))
			iSpace.(CurDataSetString){echo} = circshift(iSpace.(CurDataSetString){echo}, [0 iSpaceShift(1) 0 0 0 0 0]);
		end
	end
	if(isfield(PreProcessingInfo.(CurDataSetString), 'iSpaceShift') && numel(PreProcessingInfo.(CurDataSetString).iSpaceShift) > 1 && ne(PreProcessingInfo.(CurDataSetString).iSpaceShift(2),0))
		fprintf('\nSpatial shift data by [%s].', sprintf('%d ',[0 0 iSpaceShift(2) 0 0 0 0]))    	
		for echo = 1:numel(kSpace.(CurDataSetString))
			iSpace.(CurDataSetString){echo} = circshift(iSpace.(CurDataSetString){echo}, [0 0 iSpaceShift(2) 0 0 0 0]);
		end
	end


	if(isfield(PreProcessingInfo.(CurDataSetString), 'NoFFT_flag') && PreProcessingInfo.(CurDataSetString).NoFFT_flag)		% Remove it again
		iSpace = rmfield(iSpace,CurDataSetString);
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

memused_after = memused_linux(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


