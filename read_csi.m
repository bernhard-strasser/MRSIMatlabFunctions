function [iSpace,Noise,kSpace,PreProcessingInfo,ReadInInfo] = read_csi(file,PreProcessingInfo)
%
% read_csi Read in csi-data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function only decides if input file is .dat file or DICOM according to its ending. It then calls the read_csi_dat_x_x or read_csi_dicom_x_x
% functions. Refer to these for more info
%
%
% [csi,NoiseCorrMat,Noise_mat,csi_kspace] = read_csi(file,zerofill_to_nextpow2_flag,zerofilling_fact, Hadamard_flag, x_shift,y_shift,NoFFT_flag,NoiseCorrMat)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         zerofill_to_nextpow2_flag   ...     Flag, if the MRSI data should be zerofilled to the next power of 2 in k-space (e.g. 42x42 sampled --> zf to 64x64?)
% -         zerofilling_fact            ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% -         Hadamard_flag               ...     If data is multislice hadamard encoded, perform hadamard-decoding function
% -         x_shift                     ...     Shift the MRSI data in the left-right direction ( = row direction of matrix) by x_shift voxels
% -         y_shift                     ...     Shift the MRSI data in anterior-posterior direction ( = column direction of matrix) by y_shift voxels
% -         NoFFT_flag                  ...     If true, no fft is performed when reading dat-files.
% -         NoiseCorrMat                ...     If size(NoiseCorrMat) = [cha cha]: the k-space Data gets decorrelated with this matrix. 
%                                               If NoiseCorrMat = 1: the end of the FIDs of the outermost k-space/image points are used for noise decorrelation.
%                                               If NoiseCorrMat = 0, or not existant: No Noise Decorrelation performed
%
% Output:
% -         csi                         ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         NoiseCorrMat                ...     The Noise Correlation Matrix in order to check if it was properly computed from the csi data. Is 0 if no decorrelation performed.
% -         Noise_mat                   ...     The Noise gathered from the CSI data. Noise_mat = 0, if not gathered.
% -         csi_kspace                  ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_csi_dat_1_10, Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m, read_csi_dicom_1_1



%% 0. PREPARATIONS

% % Assign standard values to variables if nothing is passed to function.

% If nothing is passed to function
if(nargin < 1)
    display('Please feed me with a file to read in.')
    iSpace = 0;
    kSpace = 0;
    return;
end

% Test if any kSpace Preprocessing should be done with ONLINE
if(exist('PreProcessingInfo','var') && isfield(PreProcessingInfo,'ONLINE'))
	ONLINEkSpaceNecessary = ne(PreProcessingInfo.ONLINE.fredir_shift,0) && PreProcessingInfo.ONLINE.Hamming_flag;
else
	ONLINEkSpaceNecessary = false;
end



%% 1. Read In Data


if(numel(strfind(file, '.dat')) > 0)
    
    % Read Raw Data
    [kSpace,ReadInInfo] = read_csi_dat(file);
	
	
      
    
else   
   
    if(nargout == 3 || ONLINEkSpaceNecessary)
        [iSpace.ONLINE,kSpace.ONLINE] = read_csi_dicom(file);
    else
        kSpace.ONLINE = read_csi_dicom(file);      	% This is in fact iSpace, but to make it work with PreProcessMRIData_Wrapper it is just renamed wrongly
		PreProcessingInfo.ONLINE.NoFFT_flag = true;
    end
    Noise.CorrMat = 0;
    Noise.Mat = 0;
	ReadInInfo = 0;
    
end







%% 2. PreProcesing Preps


% Define Standard PreProcessingInfo
PreProcessingInfo_Standard.ONLINE.NoiseCorrMat = 1;
PreProcessingInfo_Standard.PATREFANDIMASCAN.NoiseCorrMat = 1;
PreProcessingInfo_Standard.ONLINE.Hadamard_flag = true;
PreProcessingInfo_Standard.PATREFANDIMASCAN.Hadamard_flag = true;
PreProcessingInfo_Standard.ONLINE.NoFFT_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.NoFFT_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.FlipkSpaceAlong = 2;
PreProcessingInfo_Standard.PATREFANDIMASCAN.FlipkSpaceWhileAccessing = ':,:,:,:,:,:,2';
PreProcessingInfo_Standard.ONLINE.fredir_shift = 0;
PreProcessingInfo_Standard.PATREFANDIMASCAN.fredir_shift = 0;
PreProcessingInfo_Standard.ONLINE.SaveUnfilteredkSpace = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.SaveUnfilteredkSpace = true;
PreProcessingInfo_Standard.ONLINE.Hamming_flag = true;
PreProcessingInfo_Standard.PATREFANDIMASCAN.Hamming_flag = true;


if(exist('ReadInInfo','var') && isfield(ReadInInfo.ONLINE, 'nReadEnc') && ~(ReadInInfo.ONLINE.nReadEnc == 1) )
	PreProcessingInfo_Standard.PATREFANDIMASCAN.EllipticalFilterSize = ReadInInfo.ONLINE.nReadEnc/2;
end


% If PreProcessingInfo is not passed over
if(nargin < 2)
    PreProcessingInfo = PreProcessingInfo_Standard;
end

% If PreProcessingInfo does not contain all necessary fields
PreProcessingFields = transpose(fields(kSpace)); PreProcessingFields(strcmpi(PreProcessingFields,'NOISEADJSCAN')) = [];
for i = PreProcessingFields
	if(~isfield(PreProcessingInfo,i{:}))
		PreProcessingInfo.(i{:}) = PreProcessingInfo_Standard.(i{:});
	end

	for j = transpose(fields(PreProcessingInfo_Standard.(i{:})))
		if(~isfield(PreProcessingInfo.(i{:}),j{:}))
			PreProcessingInfo.(i{:}).(j{:}) = PreProcessingInfo_Standard.(i{:}).(j{:});
		end
	end
end
clear PreProcessingFields i j;

clear PreProcessingInfo_Standard;







%% 3. PreProcess Data


% PreProcess Data
if(nargout >= 4)
	[iSpace,Noise,kSpace,PreProcessingInfo] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo);
elseif(nargout == 3)
	[iSpace,Noise,kSpace] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo);
elseif(nargout == 2)
	[iSpace,Noise] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo);
	clear kSpace;
else
	iSpace = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo);
	clear kSpace;
end






%% 3. Postparations

% if(isfield(kSpace,'Noise'))
%     if(PreProcessingInfo.Values.NoiseCorrMat == 1)
%         Noise.Mat = kSpace.Noise;
%     end
%     kSpace = rmfield(kSpace,'Noise');
% end
% if(isfield(kSpace,'NoiseCorrMat'))
%     if(PreProcessingInfo.Values.NoiseCorrMat == 1)
%         Noise.CorrMat = kSpace.NoiseCorrMat;
%     end
%     kSpace = rmfield(kSpace,'NoiseCorrMat');
% end
% 
% if(isfield(PreProcessingInfo.Values,'NoiseCorrMat'))
%     if(numel(PreProcessingInfo.Values.NoiseCorrMat) > 1)
%         Noise.CorrMat = PreProcessingInfo.Values.NoiseCorrMat;
%     end
% end





