function [kSpace, Info] = read_Siemens_dat(file)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations



Info.mdhEntryNames = {'channel','LineIndex','PhaseIndex','PartIndex','SliceIndex','EchoIndex','AcqIndex_AkaAvg','RepIndex','samples','samples before echo','SegIndex','SetIndex','ida','idb','idc','idd','FreeIcePara1','FreeIcePara2','FreeIcePara3','FreeIcePara4'};



%% 1. Gather information from header



ParList = read_ascconv(file);
Info.General.Ascconv = ParList;


% Analyze mdh.
ParList = Analyze_mdh(file,1);
%Info.ONLINE.total_channel_no = ParList.total_channel_no;
Info.General.total_ADCs = ParList.General.total_ADC_meas;




%% 2. Initialize Data

filesize = dir(file); filesize = filesize.bytes;
file_fid = fopen(sprintf('%s', file),'r');
headersize = fread(file_fid,1, 'uint32');
TotChannels = Info.General.Ascconv.total_channel_no_measured;
datapoints = ((filesize - headersize) - 128*Info.General.total_ADCs)/8/TotChannels;         % 8 bytes
% datapoints = Info.General.total_ADCs*4000;


% Allocate memory
% kSpace.(CurrentMeasSet) = zeros([ParList.(CurrentMeasSet).total_channel_no, ParList.(CurrentMeasSet).TotalRawPointsPerChannel]);
% Info.mdhInfo = zeros([18 ParList.(CurrentMeasSet).NoOfADCs*ParList.(CurrentMeasSet).total_channel_no]);

kSpace = zeros([TotChannels datapoints]);
Info.mdhInfo = zeros([20 Info.General.total_ADCs]);
Info.EvalInfo = cell([1 Info.General.total_ADCs]);


%% 3. READ DATA

tic
fprintf('\nReading data\t\t\t...')
        
fseek(file_fid, headersize,'bof'); 

chak_header = zeros([1 64]);
ACQEND_flag = false;
CurADC = 0;
CurPoint = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Loop over different measurement sets	%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(~ACQEND_flag)
    
    
    CurADC = CurADC + 1;
    
    % Read mdh
    chak_header(1:5) = fread(file_fid, 5, 'uint32');
    EvalInfoMask = fread(file_fid, 64, 'ubit1');
	CurChak = fread(file_fid, 64-14, 'int16');
    %fseek(file_fid, -128,'cof');
    
    % Stop if ACQEND was found
	if(EvalInfoMask(1) || isempty(CurChak))
        ACQEND_flag = true;
        kSpace(:,CurPoint+1:end) = [];
        continue;
	end
    
	chak_header(8:64-7) = CurChak;

	
    % Set CurrentMeasSet
	CurrentMeasSet = Associate_EvalInfoMask(EvalInfoMask); 
    

    % Save mdh-Info for later reshaping
% 		CurADC = CurADC + 1;
%     chak_header(1:18)
    Info.EvalInfo{CurADC} = CurrentMeasSet;
    Info.mdhInfo(1,CurADC) = TotChannels;                                               % channel, not useful, cause only read mdh of 1st cha
    Info.mdhInfo(2,CurADC) = chak_header(10) + 1;                                       % LineIndex
    Info.mdhInfo(3,CurADC) = chak_header(15) + 1;                                       % PhaseIndex
    Info.mdhInfo(4,CurADC) = chak_header(13) + 1;                                       % PartIndex
    Info.mdhInfo(5,CurADC) = chak_header(12) + 1;										% SliceIndex (also hada step)
    Info.mdhInfo(6,CurADC) = chak_header(14) + 1;										% EchoIndex
    Info.mdhInfo(7,CurADC) = chak_header(11) + 1;										% AcqIndex_AkaAvg
    Info.mdhInfo(8,CurADC) = chak_header(16) + 1;										% RepIndex
    Info.mdhInfo(9,CurADC) = chak_header(8);											% samples
    Info.mdhInfo(10,CurADC) = chak_header(38) *2;										% samples before echo
    Info.mdhInfo(11,CurADC) = chak_header(18) + 1;                                      % SegIndex    
    Info.mdhInfo(12,CurADC) = chak_header(17) + 1;                                      % Set-Index    
    Info.mdhInfo(13,CurADC) = chak_header(19) + 1;										% ida
    Info.mdhInfo(14,CurADC) = chak_header(20) + 1;										% idb
    Info.mdhInfo(15,CurADC) = chak_header(21) + 1;										% idc
    Info.mdhInfo(16,CurADC) = chak_header(22) + 1;										% idd
    Info.mdhInfo(17,CurADC) = chak_header(34);											% FreeIcePara1
    Info.mdhInfo(18,CurADC) = chak_header(35);											% FreeIcePara2	
    Info.mdhInfo(19,CurADC) = chak_header(36);											% FreeIcePara3
    Info.mdhInfo(20,CurADC) = chak_header(37);											% FreeIcePara4	

    
	for CurCha = 1:TotChannels
		% Read & Assign Data	% Read real & imaginary (--> Info.(CurrentMeasSet).Samples*2) measured points
		chak_data = fread(file_fid, Info.mdhInfo(9,CurADC)*2, 'float32'); 
		chak_data = chak_data(Info.mdhInfo(10,CurADC)*2+1:end); 
		chak_data = complex(chak_data(1:2:end),chak_data(2:2:end));
		kSpace(CurCha,CurPoint+1 : CurPoint+numel(chak_data)) = chak_data;
		fseek(file_fid, 128,'cof'); %Skip the mdh of the next channel (dont need mdh for each channel!)

	end
	fseek(file_fid, -128,'cof'); % Jump back to beginning of mdh, so that we can update the mdhInfo after all channels are read in
	

% %     Check if this was the last measurement of scan
% %     Temporarily disabled because the flag is set for all HadamardSteps. Has to be changed.
%     if(EvalInfoMask_loop(9))		% LASTSCANINMEAS
%         break;
%     end		


    % if(chak_header(56)+1 == Info.(CurrentMeasSet).total_channel_no)
    CurPoint = CurPoint + numel(chak_data);
    % end

end

    
fclose(file_fid);

fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)       





%% 3. Reshape Data




%% 7. Postparations



