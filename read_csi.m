function [iSpace,Noise,PreProcessingInfo,ReadInInfo,kSpace] = read_csi(file,PreProcessingInfo,ReadInDataSets)
%
% read_csi Read in csi-data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function sets default values for PreProcessingInfo and decides if input file is .dat file or DICOM according to its ending. 
% It then calls the read_csi_dat or read_csi_dicom functions. Refer to these for more info.
%
%
% [iSpace,Noise,PreProcessingInfo,ReadInInfo,kSpace] = read_csi(file,PreProcessingInfo,ReadInDataSets)
%
% Input: 
% -         file                        ...     Path of MR(S)(I) file.
% -         PreProcessingInfo           ...     Info about how the read in data should be pre-processed. See help PreProcessMRIData_Wrapper.
% -         ReadInDataSets              ...     Only read in those datasets, e.g. 'ONLINE', 'NOISEADJSCAN', 'PATREFANDIMASCAN'. These are compared to
%                                               the EvalInfoMask of the raw-data, and only if it matches, that data set is read in.
%
% Output:
% -         iSpace                      ...     Structure of data in image domain. Its fieldnames are according to the EvalInfoMask information (e.g. 'ONLINE', 'PATREFANDIMASCAN', etc.)
%                                               size of each field: channel x nFreqEnc x nPhasEnc x nPartEnc x nSlc x SpectroVecSize x Averages.
% -         Noise                       ...     Structure of Noise gathered from the Prescan or CSI data. Noise = 0, if not gathered. It is scaled individually for each measurement set.
% -         PreProcessingInfo           ...     The updated PreProcessingInfo, as it is changed when reading in.
% -         ReadInInfo                  ...     The info of the mdh.
% -         kSpace                      ...     Structure of data in k-space, see iSpace.
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

if(exist('PreProcessingInfo','var') && numel(PreProcessingInfo) == 1 && ~isstruct(PreProcessingInfo) && PreProcessingInfo == 0)
	clear PreProcessingInfo;
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
    
    % Find out Version. This code is copied from Philipp Ehses from his "mapVBVD_20150918" functions.
    file_fid = fopen(sprintf('%s', file),'r');
    firstInt  = fread(file_fid,1,'uint32');
    secondInt = fread(file_fid,1,'uint32');
    fclose(file_fid);
    if (and(firstInt < 10000, secondInt <= 64))                             % lazy software version check (VB or VD?)
        [kSpace,ReadInInfo] = read_csi_dat_VE11(file,0,ReadInDataSets);     % VE Version
    else
        [kSpace,ReadInInfo] = read_csi_dat(file,0,ReadInDataSets);          % VB Version
    end
    clear file_fid firstInt secondInt
    
else   
   
    if(nargout == 5 || ONLINEkSpaceNecessary)
        [iSpace.ONLINE{1},kSpace.ONLINE{1}] = read_csi_dicom(file);
    else
        if(exist('PreProcessingInfo','var') && isfield(PreProcessingInfo.ONLINE,'NoiseCorrMat') && numel(PreProcessingInfo.ONLINE.NoiseCorrMat) > 1)
            kSpace.ONLINE{1} = read_csi_dicom(file);      	% This is in fact iSpace, but to make it work with PreProcessMRIData_Wrapper it is just renamed wrongly. 
            PreProcessingInfo.ONLINE.NoFFT_flag = true;
        else
            iSpace.ONLINE{1} = read_csi_dicom(file);
            Noise.CorrMat = 0;
            Noise.Mat = 0;
            ReadInInfo = 0;
            return;
        end
    end
    Noise.CorrMat = 0;
    Noise.Mat = 0;
    ReadInInfo = 0;
    
end




%% 2. Define Standard PreProcessingInfo

PreProcessingInfo_Standard.ONLINE.NoiseCorrMat = 1;
PreProcessingInfo_Standard.PATREFANDIMASCAN.NoiseCorrMat = 1;
PreProcessingInfo_Standard.ONLINE.Hadamard_flag = true;				% If this flag is set to true, hadamard decoding is performed, if several slices were measure. If set to false, 
PreProcessingInfo_Standard.PATREFANDIMASCAN.Hadamard_flag = true;	% Hadamard decoding is never performed.
PreProcessingInfo_Standard.ONLINE.NoFFT_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.NoFFT_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.FlipkSpaceAlong = 2;
PreProcessingInfo_Standard.PATREFANDIMASCAN.FlipkSpaceWhileAccessing = ':,:,:,:,:,:,2';
%PreProcessingInfo_Standard.ONLINE.fredir_shift = 0;
PreProcessingInfo_Standard.PATREFANDIMASCAN.fredir_shift = 0;
PreProcessingInfo_Standard.ONLINE.SaveUnfilteredkSpace = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.SaveUnfilteredkSpace = true;
PreProcessingInfo_Standard.ONLINE.Hamming_flag = false;
PreProcessingInfo_Standard.PATREFANDIMASCAN.Hamming_flag = false;

% Remove Oversampling
if(isfield(kSpace,'ONLINE') && size(kSpace.ONLINE{1},6) > 1)
	PreProcessingInfo_Standard.ONLINE.RmOs = false;		% CSI
else
	PreProcessingInfo_Standard.ONLINE.RmOs = true;		% Imaging
end
PreProcessingInfo_Standard.PATREFANDIMASCAN.RmOs = true;
OversamplingFactor_PATREFANDIMASCAN = 2;
OversamplingFactor_ONLINE = 2;
if(isfield(ReadInInfo,'PATREFANDIMASCAN') && isfield(ReadInInfo.PATREFANDIMASCAN, 'nReadEnc') && size(kSpace.PATREFANDIMASCAN{1},2) == size(kSpace.PATREFANDIMASCAN{1},3))  % This is really really bad...
	PreProcessingInfo_Standard.PATREFANDIMASCAN.RmOs = false;
	OversamplingFactor_PATREFANDIMASCAN = 1;
end	
if(isfield(kSpace,'ONLINE') && isfield(ReadInInfo,'ONLINE') && isfield(ReadInInfo.ONLINE, 'nReadEnc') && size(kSpace.ONLINE{1},2) == size(kSpace.ONLINE{1},3))  % This is really really bad...
	PreProcessingInfo_Standard.ONLINE.RmOs = false;
	OversamplingFactor_ONLINE = 1;
end	

% Zerofilling
if(isfield(kSpace,'PATREFANDIMASCAN') && isfield(kSpace,'ONLINE') && size(kSpace.ONLINE{1},2) > 1)
	for echo = 1:numel(kSpace.ONLINE)
		bla = size(kSpace.ONLINE{echo}); bla = [bla(1:4) size(kSpace.PATREFANDIMASCAN{echo},5) 1 size(kSpace.PATREFANDIMASCAN{echo},7)];
		PreProcessingInfo_Standard.PATREFANDIMASCAN.ZeroFillingDesiredSize{echo} = bla; clear bla;
	end
end

% Elliptical Filtering
if(exist('ReadInInfo','var') && isfield(ReadInInfo, 'ONLINE') && isfield(ReadInInfo.ONLINE, 'nReadEnc') && ~(ReadInInfo.ONLINE.nReadEnc == 1) ...
	&& isfield(kSpace,'ONLINE') && size(kSpace.ONLINE{1},2) > 1)
	PreProcessingInfo_Standard.PATREFANDIMASCAN.EllipticalFilterSize = ReadInInfo.ONLINE.nReadEnc/2;
end

if(isfield(kSpace,'PATREFSCAN'))            % Concept Prescan. Cheat it to be the same as ONLINE
    PreProcessingInfo_Standard.PATREFSCAN = PreProcessingInfo_Standard.ONLINE;          
end



%% 3. Assign standard values

if(~exist('PreProcessingInfo','var'))
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




%% 4. Correct Values
if(isfield(PreProcessingInfo,'PATREFANDIMASCAN'))
	if(isfield(PreProcessingInfo.PATREFANDIMASCAN, 'EllipticalFilterSize'))
		if(numel(PreProcessingInfo.PATREFANDIMASCAN.EllipticalFilterSize) == 1)
			PreProcessingInfo.PATREFANDIMASCAN.EllipticalFilterSize = [OversamplingFactor_PATREFANDIMASCAN 1 1 PreProcessingInfo.PATREFANDIMASCAN.EllipticalFilterSize];
		else
			PreProcessingInfo.PATREFANDIMASCAN.EllipticalFilterSize(1) = OversamplingFactor_PATREFANDIMASCAN*PreProcessingInfo.PATREFANDIMASCAN.EllipticalFilterSize(1);			
		end
	end
	if(isfield(PreProcessingInfo.PATREFANDIMASCAN, 'ZeroFillingDesiredSize'))
		for echo = 1:numel(PreProcessingInfo.PATREFANDIMASCAN.ZeroFillingDesiredSize)
			PreProcessingInfo.PATREFANDIMASCAN.ZeroFillingDesiredSize{echo}(2) = PreProcessingInfo.PATREFANDIMASCAN.ZeroFillingDesiredSize{echo}(2) * OversamplingFactor_PATREFANDIMASCAN;
		end
	end
end

if(isfield(PreProcessingInfo,'ONLINE'))
	if(isfield(PreProcessingInfo.ONLINE, 'EllipticalFilterSize'))
		if(numel(PreProcessingInfo.ONLINE.EllipticalFilterSize) == 1)
			PreProcessingInfo.ONLINE.EllipticalFilterSize = [OversamplingFactor_ONLINE 1 1 PreProcessingInfo.ONLINE.EllipticalFilterSize];
		else
			PreProcessingInfo.ONLINE.EllipticalFilterSize(1) = OversamplingFactor_ONLINE*PreProcessingInfo.ONLINE.EllipticalFilterSize(1);			
		end
	end
	if(isfield(PreProcessingInfo.ONLINE, 'ZeroFillingDesiredSize'))
		PreProcessingInfo.ONLINE.ZeroFillingDesiredSize(2) = PreProcessingInfo.ONLINE.ZeroFillingDesiredSize(2) * OversamplingFactor_ONLINE;
	end
end


if(isfield(kSpace,'PATREFSCAN'))            % Concept Prescan. Cheat it to be the same as ONLINE
    PreProcessingInfo.PATREFSCAN = PreProcessingInfo.ONLINE;
end

% if(PreProcessingInfo.ONLINE.NoFFT_flag)
% 	
% 	
% 	
% end



%% 5. PreProcess Data

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






%% 6. Postparations

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





