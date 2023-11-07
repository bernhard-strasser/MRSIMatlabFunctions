function version = io_DetectSiemensTwixSoftwareVersion(file)
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


        
file_fid = fopen(file,'r');

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
else
    version = 'vb';
end

fclose(file_fid);


