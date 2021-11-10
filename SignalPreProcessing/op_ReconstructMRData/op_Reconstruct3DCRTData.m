function [MRStruct] = op_Reconstruct3DCRTData(MRStruct,Settings)
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

% Further remarks: This function uses FFTs to get from k- to MRStruct-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct();
end
if(~isfield_recursive(Settings,'ReadInTraj.GradDelay_x_us'))
    Settings.ReadInTraj.GradDelay_x_us = 0;
end
if(~isfield_recursive(Settings,'ReadInTraj.GradDelay_y_us'))
    Settings.ReadInTraj.GradDelay_y_us = 0;    
end
if(~isfield_recursive(Settings,'CalcOutTraj.fov_overgrid'))
    Settings.CalcOutTraj.fov_overgrid = 1;    
end
if(~isfield_recursive(Settings,'CalcOutTraj.OverwriteDataSize_woOvergrid'))
    Settings.CalcOutTraj.OverwriteDataSize_woOvergrid = [];    
end
if(~isfield_recursive(Settings,'NonCartReco.CircularSFTFoV_flag'))
    Settings.NonCartReco.CircularSFTFoV_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp_flag'))
    Settings.NonCartReco.DensComp_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp.AutoScale_flag'))
    Settings.NonCartReco.DensComp.AutoScale_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp.Method'))
    Settings.NonCartReco.DensComp.Method = 'ConcentricRingTrajectory_Theoretical';    
end
if(~isfield_recursive(Settings,'NonCartReco.Correct4SpatialB0_flag'))
    Settings.NonCartReco.Correct4SpatialB0_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.ConjIniSpace_flag'))
    Settings.NonCartReco.ConjIniSpace_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.ConjInkSpace_flag'))
    Settings.NonCartReco.ConjInkSpace_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.Phaseroll_flag'))
    Settings.NonCartReco.Phaseroll_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.FlipDim_flag'))
    Settings.NonCartReco.FlipDim_flag = true;    
end





%% Calc & Read Trajectories



[MRStruct] = sim_CalcCartTraj(MRStruct,Settings.CalcOutTraj);
Settings.ReadInTraj.maxR = MRStruct.OutTraj.maxR;  % Normalize spiral trajectory to maximum of Cartesian trajectory
[MRStruct] = io_ReadCRTTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);
MRStruct.Par.TargetDataSize = [size_MultiDims(MRStruct.OutTraj.GM,[3 4]) MRStruct.Par.nPartEnc*MRStruct.Par.nSLC ...
                   MRStruct.Par.vecSize MRStruct.Par.total_channel_no_measured];


% Reco MRStruct
MRStruct = op_ReconstructNonCartMRData(MRStruct,[],Settings.NonCartReco);




%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);

