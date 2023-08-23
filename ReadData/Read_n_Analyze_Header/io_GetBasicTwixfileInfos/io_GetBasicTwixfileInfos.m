function BasicInfo = io_GetBasicTwixfileInfos(MRStructOrFile)
%
% io_GetBasicTwixfileInfos Read in raw data from Siemens
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
% -         BasicInfo                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         MRStruct.mdhInfo                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.


%%
if(isstruct(MRStructOrFile) && isfield(MRStructOrFile,'DataFile'))
    DataFile = MRStructOrFile.DataFile;
else
    DataFile = MRStructOrFile;
end



%% Open File
file_fid = fopen(sprintf('%s', DataFile),'r');


% get file size
fseek(file_fid,0,'eof');
BasicInfo.fileSize = ftell(file_fid);

% start of actual measurement data (sans header)
fseek(file_fid,0,'bof');


%% Software Check and Lengths and Offsets
firstInt  = fread(file_fid,1,'uint32');
secondInt = fread(file_fid,1,'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    BasicInfo.version = 'vd';

    % number of different scans in file stored in 2nd in
    BasicInfo.NScans = secondInt;
    BasicInfo.measID = fread(file_fid,1,'uint32');
    BasicInfo.fileID = fread(file_fid,1,'uint32');
    BasicInfo.measOffset = cell(1, BasicInfo.NScans);
    BasicInfo.measLength = cell(1, BasicInfo.NScans);
    BasicInfo.IceParamOffSet = 48;
    for k=1:BasicInfo.NScans
        BasicInfo.measOffset{k} = fread(file_fid,1,'uint64');
        BasicInfo.measLength{k} = fread(file_fid,1,'uint64'); 
        fseek(file_fid, 152 - 16, 'cof');
    end
else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    BasicInfo.version  = 'vb';
    BasicInfo.measOffset{1} = 0;
    BasicInfo.measLength{1} = BasicInfo.fileSize;
    BasicInfo.NScans     = 1; % VB does not support multiple scans in one file
    BasicInfo.IceParamOffSet = 34;
end


%% Mdh Differrences between VD and VB
if(strcmpi(BasicInfo.version,'vd'))
    BasicInfo.MdhSizeInBytes = 192;
    BasicInfo.EvalInfoRelPosInBytes = 41;
    BasicInfo.ChannelMdhOffset = 8;
    BasicInfo.HeaderMdhOffset = 0;
    
else
    BasicInfo.MdhSizeInBytes = 128;
    BasicInfo.EvalInfoRelPosInBytes = 21;
    BasicInfo.ChannelMdhOffset = 0;
    BasicInfo.HeaderMdhOffset = 32;
end

%% Header Length
for CurSet = BasicInfo.NScans:BasicInfo.NScans    %NScans % Currently not implemented. Always read last scan

    cPos = BasicInfo.measOffset{CurSet};
    fseek(file_fid,cPos,'bof');
    BasicInfo.hdr_len = fread(file_fid, 1,'uint32');
    
end
BasicInfo.BeginningOfData = cPos + BasicInfo.hdr_len;


%% Close File
fclose(file_fid);

