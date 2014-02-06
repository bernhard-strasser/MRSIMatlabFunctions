function [iSpace,Noise, kSpace, PreProcessingInfo] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo)
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
% if(~exist('x_shift','var'))
%     x_shift = 0;
% end
% if(~exist('y_shift','var'))
%     y_shift = 0;
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

	CurDataSetString = CurDataSet{:};


	% 1.1 Compute Noise Correlation Matrix if necessary
	% Only read in data if NoiseCorrelationMatrix exists and is a 1
	if( isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 1 )


		%%%%%%% GET NOISE %%%%%%

		% Get data from Noise Prescan
		if(isfield(kSpace,'NOISE'))
			fprintf('\nNoise Decorrelating Using\t...\tNoise Prescan.')    
			Noise_mat = reshape(kSpace.NOISE,[size(kSpace.NOISE,1) numel(kSpace.NOISE)/size(kSpace.NOISE,1)]);

		% Or from end of FID of the outer kSpace
		elseif(isfield(kSpace,'CSI'))
			fprintf('\nNoise Decorrelating Using\t...\tCSI Data Itself.')
			
			% Gather "Noise" of the CSI data (end of FIDs in the outermost acquired circle of k-space points). Only working for 2D and Multislice
			nFIDendpoints = 400;
			TakeNPointsOutOfEnd = 200;
			randy = randperm(nFIDendpoints); randy = randy(1:TakeNPointsOutOfEnd); % Take 'em randomly
			
			% Take random points at end of FID
			Noise_csi = kSpace.CSI(:,:,:,:,:,end - (nFIDendpoints - 1) : end); Noise_csi = Noise_csi(:,:,:,:,:,randy);

			% Take only csi-points which are farest away from k-space center (so a circle with radius 31 k-space points)
			if(ReadInInfo.General.Ascconv.Full_ElliptWeighted_Or_Weighted_Acq == 2)
				[Elliptical_dummy,OuterkSpace_mask1] = EllipticalFilter_1_1(squeeze(kSpace.CSI(1,:,:,1,1,1,1)),[1 2],[1 1 1 size(kSpace.CSI,3)/2-1],1);
				[Elliptical_dummy,OuterkSpace_mask2] = EllipticalFilter_1_1(squeeze(kSpace.CSI(1,:,:,1,1,1,1)),[1 2],[1 1 1 size(kSpace.CSI,3)/2-2],1);
				OuterkSpace_mask = OuterkSpace_mask1 - OuterkSpace_mask2;
			else
				OuterkSpace_mask = zeros([size(kSpace.CSI,2), size(kSpace.CSI,3)]);
				OuterkSpace_mask(1:end,1) = 1; OuterkSpace_mask(1,1:end) = 1; OuterkSpace_mask(end,1:end) = 1; OuterkSpace_mask(1:end,end) = 1;
			end
			PI_mask = abs(squeeze(kSpace.CSI(1,:,:,1,1,1))); PI_mask(PI_mask > 0) = 1;
			csi_mask = OuterkSpace_mask .* PI_mask;
			csi_mask = repmat(logical(csi_mask), [1 1 size(kSpace.CSI,1)*size(kSpace.CSI,4)*TakeNPointsOutOfEnd]);
			csi_mask = reshape(csi_mask, [size(kSpace.CSI,2) size(kSpace.CSI,3) size(kSpace.CSI,1) size(kSpace.CSI,4) TakeNPointsOutOfEnd]);
			csi_mask = permute(csi_mask, [3 1 2 4 5]);

			Noise_mat = Noise_csi(csi_mask);
			Noise_mat = reshape(Noise_mat, [size(kSpace.CSI,1) numel(Noise_mat)/size(kSpace.CSI,1)]);  
		end

		% Copy NoiseCorrMat
		for CurDataSet2 = DataSetNames
			if(strcmpi(CurDataSet2{:},'NOISE'))
				continue;
			end

			% Rescale Noise to be the same like csi noise
			if( isfield(ReadInInfo,'NOISE') && isfield(ReadInInfo.NOISE,'Dwelltime') && isfield(ReadInInfo,CurDataSet2{:}) && isfield(ReadInInfo.(CurDataSet2{:}),'Dwelltime') )
				NoiseScalingFactor = sqrt(ReadInInfo.(CurDataSet2{:}).Dwelltime / ReadInInfo.NOISE.Dwelltime);
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
		fprintf('\nStarted')
		NoiseCorrMat_Chol = chol(PreProcessingInfo.(CurDataSetString).NoiseCorrMat/2,'lower');
		
		% Perform Noise Decorrelation Slice by Slice and Part by Part to avoid extensive memory usage 
		for SlcLoopy = 1:size(kSpace.(CurDataSetString),4)
			for PartLoopy = 1:size(kSpae.(CurDataSetString),5)
				
				TempData = kSpace.(CurDataSetString)(:,:,:,SlcLoopy,PartLoopy,:,:);
				size_TempData = size(TempData);
				TempData = reshape(TempData, [size(TempData,1) numel(TempData)/size(TempData),1]);
				fprintf('\nFirst Reshape')
				TempData = NoiseCorrMat_Chol \ TempData;
				fprintf('\nDivision')
				kSpace.(CurDataSetString)(:,:,:,SlcLoopy,PartLoopy,:,:) = reshape(TempData,size_TempData);
				fprintf('\nSecond Reshape')
				
			end
		end
		
				
				
				
% 		size_kSpace = size(kSpace.(CurDataSetString));
% 		kSpace.(CurDataSetString) = reshape(kSpace.(CurDataSetString), [size(kSpace.(CurDataSetString),1) numel(kSpace.(CurDataSetString))/size(kSpace.(CurDataSetString),1)]);
% 		fprintf('\nFirst Reshape')
% 		NoiseCorrMat_Chol = chol(PreProcessingInfo.(CurDataSetString).NoiseCorrMat/2,'lower');
% 		fprintf('\nchol')		
% 		kSpace.(CurDataSetString) = NoiseCorrMat_Chol \ kSpace.(CurDataSetString);
% 		fprintf('\nDivision')
% 		kSpace.(CurDataSetString) = reshape(kSpace.(CurDataSetString),size_kSpace);
% 		fprintf('\nSecond Reshape')

		
% 		kSpace.(CurDataSetString) = reshape(chol(PreProcessingInfo.(CurDataSetString).NoiseCorrMat/2,'lower') \ ...
% 		reshape(kSpace.(CurDataSetString), [size(kSpace.(CurDataSetString),1) numel(kSpace.(CurDataSetString))/size(kSpace.(CurDataSetString),1)]), size(kSpace.(CurDataSetString)));    % Matrix multiplication
		clear nFIDendpoints TakeNPointsOutOfEnd randy Noise_csi Elliptical_dummy OuterkSpace_mask2 OuterkSpace_mask1 OuterkSpace_mask PI_mask csi_mask
	end



	
	% 1.3. Flip
	if(isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceAlong') && isfield(PreProcessingInfo.(CurDataSetString),'FlipkSpaceWhileAccessing'))
		for FlipAlong = PreProcessingInfo.(CurDataSetString).FlipkSpaceAlong
			eval(['kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') = flipdim(kSpace.(CurDataSetString)(' PreProcessingInfo.(CurDataSetString).FlipkSpaceWhileAccessing ') ,FlipAlong);']);
			%%image_kspace = circshift(image_kspace, [0 1 1 0 0]);
		end
	end





	% 1.2 Hadamard Decoding
	if(isfield(PreProcessingInfo.(CurDataSetString),'Hadamard_flag') && PreProcessingInfo.(CurDataSetString).Hadamard_flag && size(kSpace.(CurDataSetString),5) > 1)
		for channel_no = 1:size(kSpace.CSI,1)

			fprintf('\nHadamard decoding channel %02d\t...', channel_no)
			tic

			kSpace.(CurDataSetString)(channel_no,:,:,:,:,:,:) = hadamard_decoding(kSpace.(CurDataSetString)(channel_no,:,:,:,:,:,:),5);            

			fprintf('\ttook\t%10.6f seconds', toc)              

		end
	end





	% 1.3 FFT FROM K_SPACE TO DIRECT SPACE

	if(isfield(PreProcessingInfo.(CurDataSetString), 'NoFFT_flag') && ~PreProcessingInfo.(CurDataSetString).NoFFT_flag)
		if(size(kSpace.(CurDataSetString),2) > 1 || size(kSpace.(CurDataSetString),3) > 1 || size(kSpace.(CurDataSetString),4) > 1)                    % In case that MRSI data, not Single Voxel Spectroscopy data is inputted.

			iSpace.(CurDataSetString) = kSpace.(CurDataSetString);

			if(nargout < 3)
				kSpace = rmfield(kSpace,CurDataSetString);
			end    


			% Sacher: Doing everythin at once, because enough memory is available.
			[bla, hostname] = unix('echo $(hostname)'); clear bla;
			if(~isempty(strfind(hostname,'sacher')))

					tic_overall = tic;
					iSpace.(CurDataSetString) = ifftshift(ifftshift(ifftshift(iSpace.(CurDataSetString),2),3),4);
					iSpace.(CurDataSetString) = fft(fft(fft(iSpace.(CurDataSetString),[],2),[],3),[],4);
					
					if(strcmpi(CurDataSetString,'CSI'))
						iSpace.(CurDataSetString) = conj(iSpace.(CurDataSetString));								 % the chem shift gets higher from right to left --> conj reverses that
					end
					iSpace.(CurDataSetString) = fftshift(fftshift(fftshift(iSpace.(CurDataSetString),2),3),4); 
        
					iSpace.(CurDataSetString) = flipdim(iSpace.(CurDataSetString),2);							% THIS FLIPS LEFT AND RIGHT IN SPATIAL DOMAIN BECAUSE PHYSICIANS WANT TO SEE IMAGES FLIPPED 
					fprintf('\nOverall FFT Process\t\t...\ttook\t%10.6f seconds', toc(tic_overall))

			% All others: Everything channel by channel.
			else

				tic_overall = tic;
				for channel_no = 1: size(iSpace.(CurDataSetString),1)

					tic_loop = tic;
					fprintf('\nFouriertransforming channel %02d\t...', channel_no)




					data_channel = iSpace.(CurDataSetString)(channel_no,:,:,:,:,:);                    % This extra assignment proved to be faster than using always iSpace.(CurDataSetString)(channel_no,:,:,:,:,:). Seems that indexing is rather slow.

					data_channel = ifftshift(ifftshift(ifftshift(data_channel,2),3),4);
					data_channel = fft(fft(fft(data_channel,[],2),[],3),[],4);



					data_channel = conj(data_channel);                            % the chem shift gets higher from right to left --> conj reverses that
					data_channel = fftshift(fftshift(fftshift(data_channel,2),3),4);
					data_channel = flipdim(data_channel,2);                       % THIS FLIPS LEFT AND RIGHT IN SPATIAL DOMAIN BECAUSE PHYSICIANS WANT TO SEE IMAGES FLIPPED

					iSpace.(CurDataSetString)(channel_no,:,:,:,:,:) = data_channel; 

					fprintf('\ttook\t%10.6f seconds', toc(tic_loop))       

				end

				clear data_channel
				fprintf('\nOverall FFT Process\t\t...\ttook\t%10.6f seconds', toc(tic_overall))

			end

		else
			if(strcmpi(CurDataSetString,'CSI'))
				iSpace.CSI = conj(kSpace.CSI);
				kSpace.CSI = 0;
			end
		end

	else

		iSpace.(CurDataSetString) = 0;

	end





% 	%% 6. Spatial Shift
% 
% 	if(ne(x_shift,0))
% 		iSpace.(CurDataSetString) = circshift(iSpace.(CurDataSetString), [0 x_shift 0 0 0]);
% 	end
% 
% 	if(ne(y_shift,0))
% 		iSpace.(CurDataSetString) = circshift(iSpace.(CurDataSetString), [0 0 y_shift 0 0]);
% 	end



end



%% 7. Postparations

memused_after = memused_linux_1_0(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


