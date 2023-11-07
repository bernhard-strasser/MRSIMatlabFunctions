function [MRStruct] = sim_CalcCartTraj(MRStruct,Settings)
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

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'fov_overgrid'))
    Settings.fov_overgrid = 1;
end

if(~isfield(MRStruct,'RecoPar'))
    if(~isfield(MRStruct,'Par'))
        error('Output must have field Par or RecoPar.')
    end
end
MRStruct.Par.fov_overgrid = Settings.fov_overgrid;
if(isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar.fov_overgrid = Settings.fov_overgrid;
    CurPar = MRStruct.RecoPar;
else
    CurPar = MRStruct.Par;    
end



if(isfield(Settings,'OverwriteDataSize_woOvergrid') && ~isempty(Settings.OverwriteDataSize_woOvergrid))
    DataSize = Settings.OverwriteDataSize_woOvergrid;    
else
    DataSize = [CurPar.nFreqEnc CurPar.nPhasEnc CurPar.nPartEnc*CurPar.nSLC ...
                                   CurPar.vecSize CurPar.total_channel_no_measured];
end
if(CurPar.fov_overgrid > 1)
    DataSize(1:2) = DataSize(1:2) * CurPar.fov_overgrid;   % For 3D Trajectories, this would be wrong! 
end                                                                                                          % Would also need to modify 3rd dim!

%% Calculate Cartesian Trajectory

FoV = CurPar.FoV_Read/1000*CurPar.fov_overgrid; % in m
DeltaGM = 10^9/(FoV*CurPar.GyroMagnRatioOverTwoPi);           % in mT/m * us

% Calculate a Grid for PhaseEncoding Steps
[bla_x, bla_y] = find(ones(DataSize(1)));

MRStruct.OutTraj.GM(1,1,:) = bla_x - floor(DataSize(1)/2) - 0.5; 
MRStruct.OutTraj.GM(2,1,:) = bla_y - floor(DataSize(1)/2) - 0.5; 

MRStruct.OutTraj.GM = reshape(MRStruct.OutTraj.GM * DeltaGM,[2 1 DataSize(1:2)]);

MRStruct.OutTraj.maxR = DeltaGM*DataSize(1)/2;


%% Normalize Trajectory

MRStruct.OutTraj.GM = MRStruct.OutTraj.GM/(MRStruct.OutTraj.maxR*2);


%% PostParations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);


