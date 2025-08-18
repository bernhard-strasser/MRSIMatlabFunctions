function [MRStruct,RefScan,Noise,AdditionalOut] = op_ReadAndReco3DCRTData(DataFile,TrajectoryFile,Settings)
%
% op_ReadAndReco3DCRTData Read and reconstruct data from ViennaCRT MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in ViennaCRT MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         DataFile          ...  The Siemens twix-data file of the CRT sequence   
% -         TrajectoryFile    ...  The trajectory file containing information about the k-space NonCart trajectory.
% -         Settings                ...  Struct with fields
%                                           Debug: The Debug-settings. Subfields:
%                                                           ShowTrajs:  If you want to plot the NonCart and Cartesian trajectories to get a feeling how you measured...
%                                           ReadInTraj:         Settings for reading the trajectory, see io_ReadCRTTraj.m.
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

if(~isfield_recursive(Settings,'io_Read3DCRTPars.IncludeRewinder_flag'))
   Settings.io_Read3DCRTPars.IncludeRewinder_flag = false;    
end


if(~isfield_recursive(Settings,'CalcOutTraj.fov_overgrid'))
   Settings.CalcOutTraj.fov_overgrid = 1;    
end
if(~isfield_recursive(Settings,'CalcOutTraj.OverwriteDataSize_woOvergrid'))
   Settings.CalcOutTraj.OverwriteDataSize_woOvergrid = [];    
end


if(~isfield(Settings,'ReadIn'))
    Settings.ReadIn = struct;
end
if(~isfield(Settings.ReadIn,'OmitDataSets'))
    Settings.ReadIn.OmitDataSets = {}; 
end
if(~iscell(Settings.ReadIn.OmitDataSets))
    Settings.ReadIn.OmitDataSets = {Settings.ReadIn.OmitDataSets};
end
if(nargout < 2)
	Settings.ReadIn.OmitDataSets{numel(Settings.ReadIn.OmitDataSets)+1} = 'PATREFSCAN';
    Settings.ReadIn.OmitDataSets{numel(Settings.ReadIn.OmitDataSets)+1} = 'NOISEADJSCAN';
end
if(nargout < 3)
	Settings.ReadIn.OmitDataSets{numel(Settings.ReadIn.OmitDataSets)+1} = 'NOISEADJSCAN';
end

if(~isfield(Settings,'PreWhitenData_flag'))
    Settings.PreWhitenData_flag = true;
end
if(~isfield_recursive(Settings,'NonCartReco.CircularSFTFoV_flag'))
   Settings.NonCartReco.CircularSFTFoV_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.Phaseroll_flag'))
   Settings.NonCartReco.Phaseroll_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.PhaserollRefScan_flag'))
   Settings.NonCartReco.PhaserollRefScan_flag = false;    
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
if(~isfield_recursive(Settings,'NonCartReco.DensComp.ApplyHammingFilter_flag'))
   Settings.NonCartReco.DensComp.ApplyHammingFilter_flag = false;    
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

if(~exist('Settings','var'))
    Settings = struct;
end


% Create MRStruct
if(isstruct(DataFile) && isfield(DataFile,'DataFile'))
    MRStruct = DataFile;
else
    MRStruct.DataFile = DataFile;
    if(exist('TrajectoryFile','var'))
        MRStruct.TrajFile = TrajectoryFile;
    else
        MRStruct.TrajFile = [];
    end
end

if(~isfield(MRStruct,'TrajFile'))
    MRStruct.TrajFile = [];
end

%% Find all Trajectory Files

if(~isempty(MRStruct.TrajFile) && ~iscell(MRStruct.TrajFile) && ~(endsWith(MRStruct.TrajFile,'.m') || endsWith(MRStruct.TrajFile,'.mat')))
    Files = dir(MRStruct.TrajFile);
    Files = {Files.name};
    Files = Files(endsWith(Files,'.m'));
    MRStruct.TrajFile = strcat(MRStruct.TrajFile,'/',Files);
    clear Files
end



%% Reshape, Average & Reshape Data

[MRStruct, RefScan, Noise, AdditionalOut] = io_ReadAverageReshape3DCRTDataOwnRead(MRStruct,Settings.ReadIn);


%% Reconstruct Data
if(isfield(MRStruct,'Data'))
    MRStruct = op_ReconstructMRData(MRStruct,Noise,Settings);
end
if(nargout > 1 && isfield(RefScan,'Data'))

    Settings.NonCartReco.Phaseroll_flag = Settings.NonCartReco.PhaserollRefScan_flag;
    RefScan = op_ReconstructMRData(RefScan,Noise,Settings);
end
if(nargout > 3 && isfield(AdditionalOut.CoilCompScan,'Data'))

    Settings.NonCartReco.Phaseroll_flag = Settings.NonCartReco.PhaserollRefScan_flag;    
    AdditionalOut.CoilCompScan = op_ReconstructMRData(AdditionalOut.CoilCompScan,Noise,Settings);
end



