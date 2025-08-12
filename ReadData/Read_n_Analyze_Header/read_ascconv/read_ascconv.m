function [ParList,ascconv] = read_ascconv(file_path)
%
% read_ascconv Read ascconv header part of DICOM and Siemens raw data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function cuts out the ascconv header part of DICOM and Siemens raw data and searches for Parameters within this header. These
%
%
% [ParList,ascconv] = read_ascconv(file_path)
%
% Input: 
% -         file_path                     ...     Path of file.
%
% Output:
% -         ParList                       ...     Structure giving all the Parameters. It contains among many others:
%
%           -- ParList.total_channel_no_measured         - Number of receive-channels that were measured            
%           -- ParList.total_channel_no_reco             - Number of receive-channels that were are in the file (DICOM)             
%           -- ParList.Dwelltimes                        - Dwelltime             
%           -- ParList.LarmorFreq                        - LarmorFrequency                
%           -- ParList.ThreeD_flag                       - flag telling you if measurement is real 3D measurement or 2D/Multislice                     
%           -- ParList.AsymmetricEcho                    - If AssymetricEcho was allowed                    
%           -- ParList.InterleavedSliceAcquisition       - If Slices were not measured consecutively like 1,2,3,4,... but e.g. 1,3,2,4         
%           -- ParList.nFreqEnc                          - Number of measured points in frequency encoding direction       
%           -- ParList.nPhasEnc                          - Number of measured points in phase encoding direction           
%           -- ParList.nPartEnc                          - Number of measured points in partition encoding direction (only for 3D)
%           ...
%           ...
%           ...
%
%
% -         ascconv                       ...     cell array of the ascconv header: ascconv{:,1} = ParameterName, ascconv{:,2} = ParameterValue 
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations






%%


fid = fopen(file_path, 'r', 'ieee-le');

% get file size
fseek(fid,0,'eof');
fileSize = ftell(fid);

% start of actual measurement data (sans header)
fseek(fid,0,'bof');

firstInt  = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');
fclose(fid);

if (and(firstInt < 10000, secondInt <= 64) && endsWith(file_path,'dat') && ~isdir(file_path))
    version = 'vd';
     [ParList,ascconv] = read_ascconv_VE11_eh(file_path);
else
     version = 'vb';
    [ParList,ascconv] = read_ascconv_VB(file_path);
   
end

