function MRStruct = io_ReadSiemensData(MRStructOrFile)
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
% [MRStruct.Data, MRStruct.mdhInfo] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         MRStruct.mdhInfo                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.





%% 0. Preparations

tic
fprintf('\n\tReading data\t\t\t\t...')


Settings = struct; 

% Find out memory used by MATLAB
memused_before = memused_linux(1); 

if(isstruct(MRStructOrFile))
    MRStruct = MRStructOrFile;
else
    MRStruct.DataFile = MRStructOrFile;
end
% MRStruct.mdhInfo.mdhEntryNames = {'channel','kx','ky','kz','slice','echo','avg','rep','samples','samples before echo','ida','idb','idc','idd','FreeIcePara1','FreeIcePara2','FreeIcePara3','FreeIcePara4'};

% ls DataFile, in case the user used wildcard to indicate the file (e.g. if only one dat file is in the folder, user might say file = '[path]/*.dat'

MRStruct.DataFile = strtrim(ls(MRStruct.DataFile));

%% 1. Gather information from header





        
file_fid = fopen(sprintf('%s', MRStruct.DataFile),'r');

% % Read last mdh
% fseek(file_fid, -1024,'eof');
% chak_header = fread(file_fid, 1024/2, 'uint16');
% ParList.General.total_ADC_meas = chak_header(5)-1;
% figure; plot(chak_header == 6453)
% fsfsdf

% get file size
fseek(file_fid,0,'eof');
fileSize = ftell(file_fid);

% start of actual measurement data (sans header)
fseek(file_fid,0,'bof');

firstInt  = fread(file_fid,1,'uint32');
secondInt = fread(file_fid,1,'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    version = 'vd';
    fprintf('\n\tDetected software version: VD');

    % number of different scans in file stored in 2nd in
    NScans = secondInt;
    measID = fread(file_fid,1,'uint32');
    fileID = fread(file_fid,1,'uint32');
    measOffset = cell(1, NScans);
    measLength = cell(1, NScans);
    IceParamOffSet = 48;
    for k=1:NScans
        measOffset{k} = fread(file_fid,1,'uint64');
        measLength{k} = fread(file_fid,1,'uint64'); 
        fseek(file_fid, 152 - 16, 'cof');
    end
else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    version  = 'vb';
    disp('Software version: VB (!?)');
    measOffset{1} = 0;
    measLength{1} = fileSize;
    NScans     = 1; % VB does not support multiple scans in one file
    IceParamOffSet = 34;
end

if(strcmpi(version,'vd'))
    MdhSizeInBytes = 192;
    EvalInfoRelPosInBytes = 41;
    ChannelMdhOffset = 8;
    HeaderMdhOffset = 0;
    
else
    MdhSizeInBytes = 128;
    EvalInfoRelPosInBytes = 21;
    ChannelMdhOffset = 0;
    HeaderMdhOffset = 32;
end

ParList = read_ascconv(MRStruct.DataFile);
MRStruct.Par = ParList;



%% 3. READ DATA

CurInfoPt = struct;
CurADC = struct;

for CurSet = NScans:NScans    %NScans % Currently not implemented. Always read last scan

    cPos = measOffset{CurSet};
    fseek(file_fid,cPos,'bof');
    hdr_len = fread(file_fid, 1,'uint32');
    MRStruct.Hdr = read_twix_hdr(file_fid);


    % jump to first mdh
    cPos = cPos + hdr_len;
    fseek( file_fid, cPos, 'bof' );




    chak_header = zeros([1 92]);
    ACQEND_flag = false;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%	Loop over different measurement sets	%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(~ACQEND_flag)

        % Read mdh

        CurChak = fread(file_fid, MdhSizeInBytes, 'uint8=>uint8');
        EvalInfoMask = de2bi_own(CurChak(EvalInfoRelPosInBytes:EvalInfoRelPosInBytes+7),8); EvalInfoMask = reshape(EvalInfoMask',1,[]);
        chak_header(1:5) = typecast(CurChak(1:20),'uint32');
        dummy = typecast(CurChak(EvalInfoRelPosInBytes+8:end),'uint16');
        chak_header(8:7+numel(dummy)) = dummy;

        % Stop if ACQEND was found
        if(EvalInfoMask(1) || isempty(CurChak))
            ACQEND_flag = true;
            continue;
        end


        % Set CurrentMeasSet
        CurrentMeasSet = Associate_EvalInfoMask(EvalInfoMask);
        
        if(~isfield(CurADC,CurrentMeasSet))
            CurADC.(CurrentMeasSet) = 1;
        else
            CurADC.(CurrentMeasSet) = CurADC.(CurrentMeasSet) + 1;
        end
        if(~isfield(CurInfoPt,CurrentMeasSet))
            CurInfoPt.(CurrentMeasSet) = 1;
        end
        
        %if(CurADC.NOISEADJSCAN == 34)
        %    ACQEND_flag = true;
        %    continue;        
        %end
        % Save mdh-Info for later reshaping
        MRStruct.mdhInfo.(CurrentMeasSet).Col(CurADC.(CurrentMeasSet)) = chak_header(8);                % samples
        MRStruct.mdhInfo.(CurrentMeasSet).Lin(CurADC.(CurrentMeasSet)) = chak_header(10) + 1;           % kx (Line Index)
        MRStruct.mdhInfo.(CurrentMeasSet).Ave(CurADC.(CurrentMeasSet)) = chak_header(11) + 1;           % Ave
        MRStruct.mdhInfo.(CurrentMeasSet).Sli(CurADC.(CurrentMeasSet)) = chak_header(12) + 1;           % slice (also hada step)
        MRStruct.mdhInfo.(CurrentMeasSet).Par(CurADC.(CurrentMeasSet)) = chak_header(13) + 1;           % kz (Part. Index)
        MRStruct.mdhInfo.(CurrentMeasSet).Eco(CurADC.(CurrentMeasSet)) = chak_header(14) + 1;           % echo
        MRStruct.mdhInfo.(CurrentMeasSet).Phs(CurADC.(CurrentMeasSet)) = chak_header(15) + 1;           % ky (Phase Index)
        MRStruct.mdhInfo.(CurrentMeasSet).Rep(CurADC.(CurrentMeasSet)) = chak_header(16) + 1;           % rep
        MRStruct.mdhInfo.(CurrentMeasSet).Set(CurADC.(CurrentMeasSet)) = chak_header(17) + 1;           % Set
        MRStruct.mdhInfo.(CurrentMeasSet).Seg(CurADC.(CurrentMeasSet)) = chak_header(18) + 1;           % Segment
        MRStruct.mdhInfo.(CurrentMeasSet).Ida(CurADC.(CurrentMeasSet)) = chak_header(19) + 1;           % ida
        MRStruct.mdhInfo.(CurrentMeasSet).Idb(CurADC.(CurrentMeasSet)) = chak_header(20) + 1;           % idb
        MRStruct.mdhInfo.(CurrentMeasSet).Idc(CurADC.(CurrentMeasSet)) = chak_header(21) + 1;           % idc
        MRStruct.mdhInfo.(CurrentMeasSet).Idd(CurADC.(CurrentMeasSet)) = chak_header(22) + 1;           % idd
        MRStruct.mdhInfo.(CurrentMeasSet).Ide(CurADC.(CurrentMeasSet)) = chak_header(23) + 1;           % ide
        MRStruct.mdhInfo.(CurrentMeasSet).iceParam(1:24,CurADC.(CurrentMeasSet)) = chak_header(IceParamOffSet:IceParamOffSet+23);  % FreeIcePara1-24 % I dont know why starting from 48, but it looks right
        MRStruct.mdhInfo.(CurrentMeasSet).cutoff(1:2,CurADC.(CurrentMeasSet)) = 0;	% chak_header(38) *2% samples before echo. TURNED OFF FOR NOW
	        
        
                
%         % Save mdh-Info for later reshaping
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(1,CurADC.(CurrentMeasSet)) = 1;										% channel
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(2,CurADC.(CurrentMeasSet)) = chak_header(10) + 1;		% kx (Line Index)
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(3,CurADC.(CurrentMeasSet)) = chak_header(15) + 1;		% ky (Phase Index)
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(4,CurADC.(CurrentMeasSet)) = chak_header(13) + 1;		% kz (Part. Index)
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(5,CurADC.(CurrentMeasSet)) = chak_header(12) + 1;		% slice (also hada step)
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(6,CurADC.(CurrentMeasSet)) = chak_header(14) + 1;		% echo
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(7,CurADC.(CurrentMeasSet)) = chak_header(17) + 1;		% avg
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(8,CurADC.(CurrentMeasSet)) = chak_header(16) + 1;		% rep
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(9,CurADC.(CurrentMeasSet)) = chak_header(8);			% samples
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(10,CurADC.(CurrentMeasSet)) = 0;	% chak_header(38) *2% samples before echo. TURNED OFF FOR NOW
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(11,CurADC.(CurrentMeasSet)) = chak_header(19) + 1;	% ida
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(12,CurADC.(CurrentMeasSet)) = chak_header(20) + 1;	% idb
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(13,CurADC.(CurrentMeasSet)) = chak_header(21) + 1;	% idc
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(14,CurADC.(CurrentMeasSet)) = chak_header(22) + 1;	% idd
%         MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(15:38,CurADC.(CurrentMeasSet)) = chak_header(48:71);  % FreeIcePara1-24 % I dont know why starting from 48, but it looks right
% 	      chak_data = fread(file_fid, (MRStruct.mdhInfo.(CurrentMeasSet).mdhInfo(9,CurADC.(CurrentMeasSet))*2+8)*ParList.total_channel_no_measured, 'float32'); 

        
        
        chak_data = fread(file_fid, (MRStruct.mdhInfo.(CurrentMeasSet).Col(CurADC.(CurrentMeasSet))*2+ChannelMdhOffset+HeaderMdhOffset)*ParList.total_channel_no_measured, 'float32=>float32'); 
        fseek(file_fid, -HeaderMdhOffset*4,'cof');
        if(round(numel(chak_data)/ParList.total_channel_no_measured) - numel(chak_data)/ParList.total_channel_no_measured ~= 0)
            break;
        end
        chak_data = reshape(chak_data,[numel(chak_data)/ParList.total_channel_no_measured ParList.total_channel_no_measured]);
        MRStruct.Data.(CurrentMeasSet){CurADC.(CurrentMeasSet)} = complex(chak_data(ChannelMdhOffset+1:2:end-HeaderMdhOffset,:),chak_data(ChannelMdhOffset+2:2:end-HeaderMdhOffset,:));
        


        % Check if this was the last measurement of scan
        % Temporarily disabled because the flag is set for all HadamardSteps. Has to be changed.
    %     if(EvalInfoMask_loop(9))		% LASTSCANINMEAS
    %         break;
    %     end		



    end

    if(numel(MRStruct.Data.ONLINE) == 1 && isnan(MRStruct.Data.ONLINE))
        MRStruct.Data = rmfield(MRStruct.Data,'ONLINE');
    end

    fclose(file_fid);

    fprintf('\n\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc)       

    if(numel(fields(MRStruct.Data)) > 1 && numel(fields(MRStruct.mdhInfo)) + 1 < numel(fields(MRStruct.Data)) )
        fprintf(['\nNo wipMemBlock.tFree and wipMemBlock.alFree[50-55] entries found. If you have several datasets\n(like Prescans) in your raw data, consider writing the info about\n' ...
        'the sizes of those scans in the wipMemBlock.tFree or wipMemBlock.alFree[50-55]!\n\n'])
    end

    % Create info about dimension sizes (e.g. NRep, NSet etc.)
    for CurDataset2 = transpose(fieldnames(MRStruct.Data))
        CurDataset = CurDataset2{1};
        Dummy = structfun(@numel,structfun(@unique,MRStruct.mdhInfo.(CurDataset),'uni',0),'uni',0);
        MRStruct.mdhInfo.(CurDataset) = catstruct(  MRStruct.mdhInfo.(CurDataset),cell2struct( struct2cell(Dummy),strcat('N',fieldnames(Dummy)) )  );
        MRStruct.mdhInfo.(CurDataset).NCol = MRStruct.mdhInfo.(CurDataset).Col(1);
        MRStruct.mdhInfo.(CurDataset).NCha = size(MRStruct.Data.(CurDataset){1},2);        
    end
    
    for CurDataset2 = transpose(fieldnames(MRStruct.Data))
        CurDataset = CurDataset2{1};
%         MRStruct.mdhInfo.(CurDataset).Seg = MRStruct.mdhInfo.(CurDataset).Seg - 1; % We take +1 above, but this index goes from -N/2 to N/2
        MRStruct.mdhInfo.(CurDataset).Seg(MRStruct.mdhInfo.(CurDataset).Seg > 1E4) = ...
        MRStruct.mdhInfo.(CurDataset).Seg(MRStruct.mdhInfo.(CurDataset).Seg > 1E4) - 65535 - 1;
        if(any(MRStruct.mdhInfo.(CurDataset).Seg <= 0))
        MRStruct.mdhInfo.(CurDataset).Seg = MRStruct.mdhInfo.(CurDataset).Seg + floor(MRStruct.mdhInfo.(CurDataset).NSeg/2);
        end
    end
    
    
    
end


% Comment to reshaping data:
% We cannot reshape data, because we cannot know which of the dimensions will have the same size. E.g. for CRT we can have different ADC points per
% ADC. Therefore, we can't just wirte our data as [nADCPoints nADCs nSet nAvg nSeg ...]
% The reshaping needs to be done by an own function, which knows about the sequence.





%% 7. Postparations

memused_after = memused_linux(1); 
display([' and used ' num2str(memused_after-memused_before) '% of the total memory.'])


MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);

end



