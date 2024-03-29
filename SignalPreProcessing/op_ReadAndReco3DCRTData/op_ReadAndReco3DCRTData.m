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

Noise = [];
if(nargout > 2)
    [MRStruct, RefScan, Noise] = io_ReadAverageReshape3DCRTDataOwnRead(MRStruct);
    
    standdev.RE=std(real(Noise.Data));
    standdev.IM=std(imag(Noise.Data));
    %standdev.RE=.1;
    %standdev.IM=.1;
    for ii=1:size(MRStruct.Data,2)
        MRStruct.NoiseData{ii} = single(standdev.RE*randn(size(MRStruct.Data{ii}))+1i*standdev.IM*randn(size(MRStruct.Data{ii})));
    end
    AdditionalOut.stddev = standdev; 
elseif(nargout > 1)
    [MRStruct, RefScan] = io_ReadAverageReshape3DCRTDataOwnRead(MRStruct);    
else
    MRStruct = io_ReadAverageReshape3DCRTDataOwnRead(MRStruct);    
end
% MRStruct = RefScan;
MRStruct.Data = fixnan(MRStruct.Data);

if(MRStruct.Par.Hamming_flag && Settings.NonCartReco.DensComp_flag)
    fprintf('\nWarning in op_ReadAndReco3DCRTData: Data were Hamming measured, but DensComp was on. Turned DensComp off in Reco.\n')
    Settings.NonCartReco.DensComp_flag = false;
end



%% Read Trajectory

[MRStruct] = sim_CalcCartTraj(MRStruct,Settings.CalcOutTraj);

Settings.ReadInTraj.maxR = MRStruct.OutTraj.maxR;  % Normalize NonCart trajectory to maximum of Cartesian trajectory
[MRStruct] = io_ReadCRTTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);


% % Use Measured Trajectory just to know how much to shift the trajectory (StartingPointAfterLaunchTrack),
% % but then actually use calculated and shifted trajectory
% [MRStruct2] = io_ReadCRTTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);wuff = MRStruct2.InTraj.StartingPointAfterLaunchTrack;MRStruct.TrajFile = [];
% [MRStruct] = io_ReadCRTTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);
% MRStruct.InTraj.StartingPointAfterLaunchTrack = wuff; 
% for ii = 1:numel(MRStruct.InTraj.GM); MRStruct.InTraj.GM{ii} = circshift(MRStruct.InTraj.GM{ii}, [0 -MRStruct.InTraj.StartingPointAfterLaunchTrack{ii}+1]); end


if(isfield(MRStruct.InTraj,'StartingPointAfterLaunchTrack'))
    Settings.NonCartReco.RemoveFirstADCPoints = MRStruct.InTraj.StartingPointAfterLaunchTrack;
end
    

if(nargout > 1 && isfield(RefScan,'Data'))
    [RefScan] = sim_CalcCartTraj(RefScan,Settings.CalcOutTraj);
    [RefScan] = io_ReadCRTTraj(RefScan,RefScan.TrajFile,Settings.ReadInTraj);
end

%% Spring cleaning

clearvars -except Settings MRStruct RefScan AdditionalOut Noise



%% DEBUG: PLOT NonCart Trajectories

if(Settings.Debug.ShowTrajs)
    figure;
    scatter(squeeze(MRStruct.OutTraj.GM(1,1,:)),squeeze(MRStruct.OutTraj.GM(2,1,:)),'b'), hold on   
    for AngIntNo = 1:MRStruct.Par.nAngInts
        scatter(squeeze(MRStruct.InTraj.GM{AngIntNo}(1,3:end)), squeeze(MRStruct.InTraj.GM{AngIntNo}(2,3:end)),30,'r')
        plot(squeeze(MRStruct.InTraj.GM{AngIntNo}(1,:)), squeeze(MRStruct.InTraj.GM{AngIntNo}(2,:)),'r')
        scatter(squeeze(MRStruct.InTraj.GM{AngIntNo}(1,1)), squeeze(MRStruct.InTraj.GM{AngIntNo}(2,1)),30,'m')
        scatter(squeeze(MRStruct.InTraj.GM{AngIntNo}(1,2)), squeeze(MRStruct.InTraj.GM{AngIntNo}(2,2)),	30,'g')

    end
    hold off
end



%% Reconstruct NonCart Data

% Remark: Currently, the CRT data is a cell, for each circle one cell. This is to save the data efficiently, because different circles
% can have different trajectory points
[MRStruct,AdditionalOut.RecoOperators] = op_ReconstructNonCartMRData(MRStruct,[],Settings.NonCartReco);
if(nargout > 1 && isfield(RefScan,'Data'))
    Bak = Settings.NonCartReco.Phaseroll_flag;
    Settings.NonCartReco.Phaseroll_flag = Settings.NonCartReco.PhaserollRefScan_flag;
    RefScan = op_ReconstructNonCartMRData(RefScan,[],Settings.NonCartReco);
    Settings.NonCartReco.Phaseroll_flag = Bak;
end

%% Slice Reco

MRStruct = op_SliceReco(MRStruct);
if(nargout > 1 && isfield(RefScan,'Data'))
    RefScan = op_SliceReco(RefScan);
end


