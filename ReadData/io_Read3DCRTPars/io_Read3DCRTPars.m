function [MRStruct] = io_Read3DCRTPars(MRStruct,Settings)
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
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked MRStruct option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations


if(~exist('MRStruct','var'))
    MRStruct = struct;
end
if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isfield(Settings,'IncludeRewinder_flag'))
    Settings.IncludeRewinder_flag = false;
end

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



%% Read Ascconv

if(strcmpi(MRStruct.MapVBVDObj.softwareVersion,'vd'))
    MRStruct.Par = read_ascconv_VE11_eh(MRStruct.DataFile);
else
    MRStruct.Par = read_ascconv(MRStruct.DataFile);
end


%% Define & Calculate Other Parameters


    
% Definitions
MRStruct.Par.nAngInts = MRStruct.MapVBVDObj.NLin;


% Calculate nTIs for all circles
MRStruct.Par.nTempIntsPerAngInt = zeros([1 MRStruct.Par.nAngInts]); 
MRStruct.Par.nPartEncsPerAngInt = zeros([MRStruct.MapVBVDObj.NSeg MRStruct.Par.nAngInts]); 
MRStruct.Par.nADCsPerAngInt = zeros([1 MRStruct.Par.nAngInts]);
MRStruct.Par.nPtsPerADC = zeros([1 MRStruct.Par.nAngInts]);
MRStruct.Par.TrajPts = zeros([1 MRStruct.Par.nAngInts]);
for CurCrcl = 1:MRStruct.Par.nAngInts
    MRStruct.Par.nTempIntsPerAngInt(CurCrcl) = max(MRStruct.MapVBVDObj.Idb(MRStruct.MapVBVDObj.Lin == CurCrcl));        % Only works for 2D
    MRStruct.Par.nPartEncsPerAngInt(CurCrcl) = unique(MRStruct.MapVBVDObj.Seg(MRStruct.MapVBVDObj.Lin == CurCrcl)) + floor(MRStruct.MapVBVDObj.NSeg/2);    % e.g. from (-3,-2,-1,0,1,2,3) --> (0,1,2,3,4,5,6,7). However, mapVBVD adds +1                                                                                            
                                                                                            % so it should be (1,2,3,4,5,6,7,8).
    MRStruct.Par.nADCsPerAngInt(CurCrcl) = max(MRStruct.MapVBVDObj.Ida(MRStruct.MapVBVDObj.Lin == CurCrcl));     
    MRStruct.Par.nPtsPerADC(CurCrcl) = max(MRStruct.MapVBVDObj.Col(MRStruct.MapVBVDObj.Lin == CurCrcl));
    MRStruct.Par.TrajPts(CurCrcl) = 2*max(MRStruct.MapVBVDObj.iceParam(6,MRStruct.MapVBVDObj.Lin == CurCrcl));     % 2 because of ADC-oversampling!

end

if(any(MRStruct.Par.TrajPts == 0))
    MRStruct.Par.TrajPts = repmat(2*(MRStruct.MapVBVDObj.Idc(1)-1),[1 MRStruct.Par.nAngInts]);
end

% Define which angular interleaves are measured for each kz (matrix of nPartEncmax x nCircles)
% MRStruct.Par.nPartEncsPerAngInt = [9 9 9 7 7 7 7 5 5 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1]; % For simulating 3D
dummy1 = max(MRStruct.Par.nPartEncsPerAngInt);
dummy2 = ceil(dummy1/2);
for kz = 1:max(MRStruct.Par.nPartEncsPerAngInt)
    MRStruct.Par.AngIntsPerPartEnc(kz,:) = MRStruct.Par.nPartEncsPerAngInt >= dummy1-(     (kz - (kz>dummy2)*2*mod(kz,dummy2))        -1)*2;
end
    
% Hack -- Delete me later
if(MRStruct.Par.TrajPts == 0)
    MRStruct.Par.TrajPts = [4 10 18 24 30 36 45 48 60 60 72 80 80 90 120 120 120 120 120 144 144 144 144 180 180 180 180 180 180 216 216 216] * 2;
end



MRStruct.Par.RewPts = 0;
MRStruct.Par.TrajTotPts = MRStruct.Par.TrajPts + MRStruct.Par.RewPts;



MRStruct.Par.VTI_Flag = ~all(MRStruct.Par.nTempIntsPerAngInt == MRStruct.Par.nTempIntsPerAngInt(1));



% For getting TrajPts, it depends on whether we read in ONLINE or refscan data
if(strcmpi(MRStruct.MapVBVDObj.dataType,'image'))
     MRStruct.Par.vecSize = MRStruct.MapVBVDhdr.Dicom.alICEProgramPara(7);     
else
    MRStruct.Par.vecSize = MRStruct.MapVBVDhdr.Dicom.alICEProgramPara(8);
    if(any(mod(MRStruct.Par.vecSize,MRStruct.Par.nTempIntsPerAngInt)))
        MRStruct.Par.vecSize = MRStruct.MapVBVDhdr.Dicom.alICEProgramPara(8) * max(MRStruct.Par.nTempIntsPerAngInt);    % In one version of my sequence we need to do that...
    end
    
end
if(MRStruct.Par.vecSize <= 0)
    % What's that for?
%     if(MRStruct.Par.VTI_Flag)
%         MRStruct.Par.vecSize = floor(MRStruct.Par.vecSize  / factorial(max(MRStruct.Par.nTempIntsPerAngInt))) * factorial(max(MRStruct.Par.nTempIntsPerAngInt)); 
%     end
    MRStruct.Par.vecSize = round(MRStruct.Par.nPtsPerADC(1)*MRStruct.Par.nADCsPerAngInt(1)/MRStruct.Par.TrajPts(1)*MRStruct.Par.nTempIntsPerAngInt(1)-0.5);
end
MRStruct.Par.UsefulADCPtsPerAngInt=MRStruct.Par.vecSize./MRStruct.Par.nTempIntsPerAngInt.*MRStruct.Par.TrajPts;
MRStruct.Par.ADCdtPerAngInt_ns = MRStruct.Par.Dwelltimes(1).*MRStruct.Par.nTempIntsPerAngInt./MRStruct.Par.TrajPts;       
    % Factor 2 NO LONGER NECESSARY, it's done in read_ascconv_VE11_eh!

MRStruct.Par.SpatialSpectralEncoding_flag = true;
MRStruct.Par.SpatialSpectralEncoding_type = 'CRT';


%% Fix Pars
MRStruct = supp_FixPars(MRStruct);


