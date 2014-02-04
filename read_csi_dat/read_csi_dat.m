function kSpace = read_csi_dat(csi_path, DesiredSize)
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
% [csi,NoiseCorrMat,Noise_mat,kSpace] = read_csi_dat(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag, NoiseCorrMat)
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
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
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



if(exist('DesiredSize','var') && ~isstruct(DesiredSize))
	DesiredSize2 = DesiredSize; clear DesiredSize; DesiredSize.CSI = DesiredSize2; clear DesiredSize2;
end






%% Gather information about sizes


ParList = read_ascconv(csi_path);
% First try: Get info from wipMemBlock_tfree
if(ParList.wipMemBlock_tFree ~= 0)
	% Only evaluate the tFree string if it is proper MATLAB CODE (doesnt throw any errors)
	try
		eval(ParList.wipMemBlock_tFree);
	catch errie %#ok
		if(exist('errie','var'))
			fprintf(['\nThe wipMemBlock.tFree entry is not MATLAB code. Consider writing the info about\n' ...
			'the sizes of those scans in the wipMemBlock.tFree!\n\n'])
			clear errie;
			ParList.wipMemBlock_tFree = 0;
		end
	end
end

% Second try: Try to get Prescans Info from wipMemBlock and CSI Info from normal ascconv header and mdh.
% Prescan Info
if(ParList.wipMemBlock_alFree ~= 0)
	try
		Info.GRE.nReadEnc = ParList.wipMemBlock_alFree(1);
		Info.GRE.nPhasEnc = ParList.wipMemBlock_alFree(2);
		Info.GRE.nPartEnc = ParList.wipMemBlock_alFree(3);
		Info.GRE.nSLC = ParList.wipMemBlock_alFree(4);
		Info.GRE.nAverages = ParList.wipMemBlock_alFree(5);
		Info.NOISE.nReadEnc = ParList.wipMemBlock_alFree(6);
	catch errie %#ok
		if(exist('errie','var'))
			fprintf(['\nIs there something different than infos about the Prescans in wipMemBlock.alFree[50-55] ?\n' ...
			'Consider writing the Prescans Info into that array for faster read-in.'])
			clear errie;
			ParList.wipMemBlock_alFree = 0;
		end
	end
end
	
% CSI Info
if(ParList.wipMemBlock_tFree == 0)
	if(ParList.wipMemBlock_alFree == 0)
		fprintf(['\nNo wipMemBlock.tFree and wipMemBlock.alFree[50-55] entries found. If you have several datasets\n(like Prescans) in your raw data, consider writing the info about\n' ...
		'the sizes of those scans in the wipMemBlock.tFree!\n\n'])
	end
	Info.CSI.nReadEnc = ParList.nFreqEnc;
	Info.CSI.nPhasEnc = ParList.nPhasEnc;
	Info.CSI.nPartEnc = ParList.nPartEnc;
	Info.CSI.nSLC = ParList.nSLC;
	Info.CSI.nAverages = ParList.nAverages;
end


% Analyze mdh.
ParList = Analyze_csi_mdh(csi_path,0);
%Info.CSI.total_channel_no = ParList.total_channel_no;
Info.General.total_ADCs = ParList.total_ADC_meas;
clear ParList;


% Otherwise read data without memory preallocation. This can be slow.




%% 2. READ DATA

tic
fprintf('\nReading data\t\t\t...')
        
csi_fid = fopen(sprintf('%s', csi_path),'r');
headersize = fread(csi_fid,1, 'uint32');
fseek(csi_fid, headersize,'bof'); 

chak_header = zeros([1 64]);
kSpace = struct('CSI',NaN);
ACQEND_flag = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Loop over different measurement sets	%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(~ACQEND_flag)

    % Read first mdh
    chak_header(1:5) = fread(csi_fid, 5, 'uint32');
    chak_header(6:7) = fread(csi_fid, 2, 'uint32');
    chak_header(8:64-7) = fread(csi_fid, 64-14, 'int16');
    fseek(csi_fid, -128,'cof');
    EvalInfoMask = Interpret_EvalInfoMask(chak_header(6:7));
    
    % Stop if ACQEND was found
    if(EvalInfoMask(1))
        ACQEND_flag = true;
        continue;
    end
    
    % Set CurrentMeasSet
	CurrentMeasSet = Associate_EvalInfoMask(EvalInfoMask);
	
	
	% Initialize Current Info
	if(~isfield(Info,CurrentMeasSet))
		Info.(CurrentMeasSet) = struct('nReadEnc',1,'nPhasEnc',1,'nPartEnc',1,'nSLC',1, 'nAverages',1);
	end
	
	dim = 0;
	for fieldly = {'nReadEnc','nPhasEnc','nPartEnc','nSLC','nAverages'}
		dim = dim + 1;
		% Initialize missing fields with ones
		if(~isfield(Info.(CurrentMeasSet), fieldly{:}))
			Info.(CurrentMeasSet).(fieldly{:}) = 1;
			DoShift = false;
		else
			DoShift = true;
		end		
		if(DoShift && exist('DesiredSize','var') && isfield(DesiredSize,CurrentMeasSet) && numel(DesiredSize.(CurrentMeasSet)) >= dim)
			kSpaceShift.(CurrentMeasSet)(dim) = floor(Info.(CurrentMeasSet).(fieldly{:})/2) - floor(DesiredSize.(CurrentMeasSet)(dim)/2);
			Info.(CurrentMeasSet).(fieldly{:}) = DesiredSize.(CurrentMeasSet)(dim);			
		else
			kSpaceShift.(CurrentMeasSet)(dim) = 0;
			DesiredSize.(CurrentMeasSet)(dim) = 99999;		% So that it has no effect			
		end
		
	end
	

	Info.(CurrentMeasSet).total_channel_no = chak_header(9);
	Info.(CurrentMeasSet).Samples = chak_header(8);
	
	if( strcmpi(CurrentMeasSet,'CSI') || strcmpi(CurrentMeasSet,'NOISE') )
		vecSize = Info.(CurrentMeasSet).Samples;
	else
		vecSize = 1;
		Info.(CurrentMeasSet).nReadEnc = Info.(CurrentMeasSet).Samples;
	end
	
	
	% Allocate memory
	if( ~isfield(kSpace,CurrentMeasSet) || isnan(kSpace.(CurrentMeasSet)) )					% By this statement, interleaved Prescan measurements are not allowed (data will be overwritten)
		kSpace.(CurrentMeasSet) = zeros([Info.(CurrentMeasSet).total_channel_no,Info.(CurrentMeasSet).nReadEnc,Info.(CurrentMeasSet).nPhasEnc, ...
		Info.(CurrentMeasSet).nPartEnc,Info.(CurrentMeasSet).nSLC,vecSize,Info.(CurrentMeasSet).nAverages]);
	end
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%	Loop over all measurements	%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BreakOutOfPrison = false;
    for ADC_MeasNo = 1 : Info.General.total_ADCs;
		
		if(BreakOutOfPrison)
			break;
		end
		
		for channel_no = 1:Info.(CurrentMeasSet).total_channel_no
			
			% Read again EvalInfoMask.
            chak_header(1:7) = fread(csi_fid, 7, 'uint32');
			
			% If the EvalInfoMask has changed within the loop, break out of loops.
			EvalInfoMask_loop = Interpret_EvalInfoMask(chak_header(6:7));
			CurrentMeasSet_loop = Associate_EvalInfoMask(EvalInfoMask_loop);
			if(~strcmpi(CurrentMeasSet,CurrentMeasSet_loop))
				fseek(csi_fid,-7*4,'cof');
				BreakOutOfPrison = true;
				break;
			end
			
			% Read rest of the mdh
			chak_header(8:64-7) = fread(csi_fid, 64-14, 'int16');
            k_x = chak_header(17-7) + 1 - kSpaceShift.(CurrentMeasSet)(1);                     
            k_y = chak_header(22-7) + 1 - kSpaceShift.(CurrentMeasSet)(2);  
			k_z = chak_header(20-7) + 1 - kSpaceShift.(CurrentMeasSet)(3);
            slice = chak_header(19-7) + 1;                                % SAYS WHICH REPETITION FOR HADAMARD ENCODING OF THE SAME K-POINT IS MEASURED
			Avg = chak_header(24-7) + 1;									% Averages
            %channel_no = chak_header(63-7) + 1;							% Problematic if Channel IDs are not consecutive (1,2,3,...,8 e.g., but 1,2,3,4,11,12,13,14)
            
			% Check if the k-points are alright
			if(k_x < 1 || k_x > DesiredSize.(CurrentMeasSet)(1))
				if(k_x < 1)
					fprintf('\nProblem: Detected mdh with kx < 1. Ignoring this ADC.')
				end
                fseek(csi_fid,Info.(CurrentMeasSet).Samples*2*4,'cof');
                continue;
			end
			if(k_y < 1 || k_y > DesiredSize.(CurrentMeasSet)(2))
				if(k_y < 1)
					fprintf('\nProblem: Detected mdh with ky < 1. Ignoring this ADC.')
				end
                fseek(csi_fid,Info.(CurrentMeasSet).Samples*2*4,'cof');                
                continue;
			end
			if(k_z < 1 || k_z > DesiredSize.(CurrentMeasSet)(3))
				if(k_z < 1)
					fprintf('\nProblem: Detected mdh with kz < 1. Ignoring this ADC.')
				end
                fseek(csi_fid,Info.(CurrentMeasSet).Samples*2*4,'cof');                
                continue;
			end
			if(slice < 1)
                fprintf('\nProblem: Detected mdh with kz < 1. Setting kz = 1.')
                slice = 1;															% some distinction between multislice/hadamard and real 3d-CSI necessary!
			end
			
			
			% Read & Assign Data
            chak_data = fread(csi_fid, Info.(CurrentMeasSet).Samples*2, 'float32'); % Read real & imaginary (--> Info.(CurrentMeasSet).Samples*2) measured points
			if( strcmpi(CurrentMeasSet,'CSI') || strcmpi(CurrentMeasSet,'NOISE') )
                kSpace.(CurrentMeasSet)(channel_no,k_x,k_y,k_z,slice,:,Avg) = complex(chak_data(1:2:end),chak_data(2:2:end));
            else
                kSpace.(CurrentMeasSet)(channel_no,:,k_x,k_z,slice,1,Avg) = complex(chak_data(1:2:end),chak_data(2:2:end));
			end
			
			% Check if this was the last measurement of scan
			% Temporarily disabled because the flag is set for all HadamardSteps. Has to be changed.
% 			if(EvalInfoMask_loop(9))		% LASTSCANINMEAS
% 				BreakOutOfPrison = true;
% 				break;
% 			end		
			

		end 
    end

end
    
    
fclose(csi_fid);

fprintf('\ttook\t%10.6f seconds',toc)       




%% 7. Postparations

memused_after = memused_linux_1_0(1); 
display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


