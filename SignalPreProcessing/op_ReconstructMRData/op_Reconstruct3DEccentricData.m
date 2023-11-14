function [MRStruct] = op_Reconstruct3DEccentricData(MRStruct,Settings)
%
% op_ReadAndReco3DEccentricData Read and reconstruct data from ViennaEccentric MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in ViennaEccentric MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         DataFile          ...  The Siemens twix-data file of the Eccentric sequence   
% -         TrajectoryFile    ...  The trajectory file containing information about the k-space NonCart trajectory.
% -         Settings                ...  Struct with fields
%                                           Debug: The Debug-settings. Subfields:
%                                                           ShowTrajs:  If you want to plot the NonCart and Cartesian trajectories to get a feeling how you measured...
%                                           ReadInTraj:         Settings for reading the trajectory, see io_ReadEccentricTraj.m.
%                                           CalcOutTraj:        Settings for calculating the Cartesian trajectory, see sim_CalcCartTraj.m.
%                                           NonCartReco:        Settings for the non-Cartesian MRSI Reco, see op_ReconstructNonCartMRData.m.
%
% Output:
% -         MRStruct            ...  Output-struct with fields
%                                           RecoSteps:  All the settings of the performed reconstruction steps
%                                           Par:        The original parameters of the read in data
%                                           RecoPar:    The parameters after reconstruction
%                                           Data:       The reconstructed data
%                                           NoiseData:  If available, noise with the same scale and size as the Data is produced. Useful for SNR calculations.
%                                           OutTraj:    The Cartesian trajectory which would correspond to the Fourier transform of the Cartesian output image
%                                                       (if we did k-space gridding, those would be the k-space points we would grid to).
%                                           InTraj:     The NonCart trajectory with which the data were measured.
%                                           
% -         AdditionalOut           ...  Struct with additional output, e.g. Fourier Transform Operator, Density Compensation Function, etc.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ???

% Further remarks:



%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end

if(~isfield_recursive(Settings,'Debug.ShowTrajs'))
    Settings.Debug.ShowTrajs = false;
end

if(~isfield_recursive(Settings,'io_Read3DEccentricPars.IncludeRewinder_flag'))
   Settings.io_Read3DEccentricPars.IncludeRewinder_flag = false;    
end


% if(~isfield_recursive(Settings,'ReadInTraj.IncludeRewinder_flag'))
%    Settings.ReadInTraj.IncludeRewinder_flag = false;    
% end

if(~isfield_recursive(Settings,'CalcOutTraj.fov_overgrid'))
   Settings.CalcOutTraj.fov_overgrid = 1;    
end
if(~isfield_recursive(Settings,'CalcOutTraj.OverwriteDataSize_woOvergrid'))
   Settings.CalcOutTraj.OverwriteDataSize_woOvergrid = [];    
end


if(~isfield_recursive(Settings,'NonCartReco.CircularSFTFoV_flag'))
   Settings.NonCartReco.CircularSFTFoV_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.Phaseroll_flag'))
   Settings.NonCartReco.Phaseroll_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp_flag'))
   Settings.NonCartReco.DensComp_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp.AutoScale_flag'))
   Settings.NonCartReco.DensComp.AutoScale_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp.Method'))
%    Settings.NonCartReco.DensComp.Method = 'AntoinesVoronoi'; 
   Settings.NonCartReco.DensComp.Method = 'Hamming';
%    Settings.NonCartReco.DensComp.Method = 'SpiralHoge1997AbruptChanges'; 
end
if(~isfield_recursive(Settings,'NonCartReco.ConjInkSpace_flag'))
   Settings.NonCartReco.ConjInkSpace_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.ConjIniSpace_flag'))
   Settings.NonCartReco.ConjIniSpace_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.Correct4SpatialB0_flag'))
   Settings.NonCartReco.Correct4SpatialB0_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.FlipDim_flag'))
    Settings.NonCartReco.FlipDim_flag = true;
end
if(~isfield_recursive(Settings,'NonCartReco.PseudoReplica_flag'))
   Settings.NonCartReco.PseudoReplica_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.nFIDpts'))
    Settings.NonCartReco.nFIDpts = 1;
end

if(~exist('Settings','var'))
    Settings = struct;
end


% % Create MRStruct
% if(isstruct(DataFile) && isfield(DataFile,'DataFile') && isfield(DataFile,'TrajectoryFile'))
%     MRStruct = DataFile;
% else
%     MRStruct.DataFile = DataFile;
%     MRStruct.TrajFile = TrajectoryFile;
% end

% Create MRStruct
if(~isfield(MRStruct,'TrajFile'))
    MRStruct.TrajFile = [];
end


%% Calc & Read Trajectories



[MRStruct] = sim_CalcCartTraj(MRStruct,Settings.CalcOutTraj);
Settings.ReadInTraj.maxR = MRStruct.OutTraj.maxR;  % Normalize spiral trajectory to maximum of Cartesian trajectory
[MRStruct] = io_ReadEccentricTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);
MRStruct.Par.TargetDataSize = [size_MultiDims(MRStruct.OutTraj.GM,[3 4]) MRStruct.Par.nPartEnc*MRStruct.Par.nSLC ...
                   MRStruct.Par.vecSize MRStruct.Par.total_channel_no_measured];


if(isfield(MRStruct.InTraj,'StartingPointAfterLaunchTrack'))
    Settings.NonCartReco.RemoveFirstADCPoints = MRStruct.InTraj.StartingPointAfterLaunchTrack;
end

% Reco MRStruct
MRStruct = op_ReconstructNonCartMRData(MRStruct,[],Settings.NonCartReco);

%% Slice Reco

MRStruct = op_SliceReco(MRStruct);

%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);


