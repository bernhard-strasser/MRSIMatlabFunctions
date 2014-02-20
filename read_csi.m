function [iSpace,Noise,PreProcessingInfo,ReadInInfo,kSpace] = read_csi(file,PreProcessingInfo,ReadInDataSets)
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
% -         file                        ...     Path of MR(S)(I) file.
% -         PreProcessingInfo           ...     Info about how the read in data should be pre-processed. See help PreProcessMRIData_Wrapper.
%
% Output:
% -         iSpace                      ...     Structure of data in image domain. Its fieldnames are according to the EvalInfoMask information (e.g. 'ONLINE', 'PATREFANDIMASCAN', etc.)
%                                               size of each field: channel x nFreqEnc x nPhasEnc x nPartEnc x nSlc x SpectroVecSize x Averages.
% -         Noise                       ...     Structure of Noise gathered from the Prescan or CSI data. Noise = 0, if not gathered. It is scaled individually for each measurement set.
% -         kSpace                      ...     Structure of data in k-space.
%                                               size of each field: channel x nFreqEnc x nPhasEnc x nPartEnc x nSlc x SpectroVecSize x Averages.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_csi_dat, read_csi_dicom, PreProcessMRIData_Wrapper, and many more (???)
%
% See also read_csi_dat, read_csi_dicom, PreProcessMRIData_Wrapper.



%% 0. PREPARATIONS

% % Assign standard values to variables if nothing is passed to function.

% If nothing is passed to function
if(nargin < 1)
    display('Please feed me with a file to read in.')
    iSpace = 0;
    kSpace = 0;
    return;
end
if(~exist('ReadInDataSets','var'))
	ReadInDataSets = 'All';
end

% Test if any kSpace Preprocessing should be done with ONLINE
if(exist('PreProcessingInfo','var') && isfield(PreProcessingInfo,'ONLINE') &&  isfield(PreProcessingInfo.ONLINE,'Hamming_flag') )
	ONLINEkSpaceNecessary = PreProcessingInfo.ONLINE.Hamming_flag;
else
	ONLINEkSpaceNecessary = false;
end



%% 1. Read In Data


if(numel(strfind(file, '.dat')) > 0)
    
    % Read Raw Data
    [kSpace,ReadInInfo] = read_csi_dat(file,0,ReadInDataSets);
	
	
      
    
else   
   
    if(nargout == 5 || ONLINEkSpaceNecessary)
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
PreProcessingInfo_Standard.ONLINE.Hadamard_flag = true;				% If this flag is set to true, hadamard decoding is performed, if several slices were measure. If set to false, 
PreProcessingInfo_Standard.PATREFANDIMASCAN.Hadamard_flag = true;	% Hadamard decoding is never performed.
PreProcessingInfo_Standard.ONLINE.NoFFT_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.NoFFT_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.FlipkSpaceAlong = 2;
PreProcessingInfo_Standard.PATREFANDIMASCAN.FlipkSpaceWhileAccessing = ':,:,:,:,:,:,2';
PreProcessingInfo_Standard.ONLINE.fredir_shift = 0;
PreProcessingInfo_Standard.PATREFANDIMASCAN.fredir_shift = 0;
PreProcessingInfo_Standard.ONLINE.SaveUnfilteredkSpace = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.SaveUnfilteredkSpace = true;
PreProcessingInfo_Standard.ONLINE.Hamming_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.Hamming_flag = false;



if(isfield(kSpace,'PATREFANDIMASCAN'))
	bla = size(kSpace.ONLINE); bla = [bla(1:5) 1 size(kSpace.PATREFANDIMASCAN,7)];
	PreProcessingInfo_Standard.PATREFANDIMASCAN.ZeroFillingDesiredSize = bla; clear bla;
end
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
clear PreProcessingFields i j PreProcessingInfo_Standard;







%% 3. PreProcess Data


% PreProcess Data
if(nargout >= 5)
	[iSpace,Noise,PreProcessingInfo,kSpace] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo);
elseif(nargout >= 3)
	[iSpace,Noise,PreProcessingInfo] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo);
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





