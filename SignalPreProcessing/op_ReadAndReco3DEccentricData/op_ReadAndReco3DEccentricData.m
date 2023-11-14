function [MRStruct,RefScan,AdditionalOut] = op_ReadAndReco3DEccentricData(DataFile,TrajectoryFile,Settings,standdev)
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
   Settings.NonCartReco.CircularSFTFoV_flag = false;    
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
   Settings.NonCartReco.DensComp.Method = 'AntoinesVoronoi'; 
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



if(~exist('standdev','var') && Settings.NonCartReco.PseudoReplica_flag)
    fprintf('\nWARNING in op_ReadAndReco3DEccentricData: PseudoReplica_flag was 1, althoguh no standdev variable was provided.\nSet PseudoReplica_flag to 0.');
    Settings.NonCartReco.PseudoReplica_flag = false;    
end
% Create MRStruct
if(isstruct(DataFile) && isfield(DataFile,'DataFile') && isfield(DataFile,'TrajectoryFile'))
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

if(nargout > 1)
    [MRStruct, RefScan] = io_ReadAverageReshape3DEccentricDataOwnRead(MRStruct);
else
MRStruct = io_ReadAverageReshape3DEccentricDataOwnRead(MRStruct);    
end


if (Settings.NonCartReco.PseudoReplica_flag)
    for ii=1:size(MRStruct.Data,2)
            temp = MRStruct.Data{ii}(:,1,1,1,1:Settings.NonCartReco.nFIDpts); 
            clear MRStruct.Data{ii}

            for jj=Settings.NonCartReco.nFIDpts+1:MRStruct.Par.vecSize
               MRStruct.Data{ii}(:,jj) = single(standdev.RE*randn(MRStruct.Par.TrajPts(1),1)+1i*standdev.IM*randn(MRStruct.Par.TrajPts(1),1));
            end
            clear temp 
    end
end
%% Read Trajectory


[MRStruct] = sim_CalcCartTraj(MRStruct,Settings.CalcOutTraj);

Settings.ReadInTraj.maxR = MRStruct.OutTraj.maxR;  % Normalize NonCart trajectory to maximum of Cartesian trajectory
[MRStruct] = io_ReadEccentricTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);

% % Use Measured Trajectory just to know how much to shift the trajectory (StartingPointAfterLaunchTrack), 
% % but then actually use calculated and shifted trajectory.
% [MRStruct2] = io_ReadEccentricTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);wuff = MRStruct2.InTraj.StartingPointAfterLaunchTrack;MRStruct.TrajFile = [];
% [MRStruct] = io_ReadEccentricTraj(MRStruct,MRStruct.TrajFile,Settings.ReadInTraj);
% MRStruct.InTraj.StartingPointAfterLaunchTrack = wuff; 
% for ii = 1:numel(MRStruct.InTraj.GM); MRStruct.InTraj.GM{ii} = circshift(MRStruct.InTraj.GM{ii}, [0 -MRStruct.InTraj.StartingPointAfterLaunchTrack{ii}+1]); end
    

if(isfield(MRStruct.InTraj,'StartingPointAfterLaunchTrack'))
    Settings.NonCartReco.RemoveFirstADCPoints = MRStruct.InTraj.StartingPointAfterLaunchTrack;
end


%% Spring cleaning

clearvars -except Settings MRStruct RefScan standdev



%% Reconstruct NonCart Data

% Remark: Currently, the Eccentric data is a cell, for each circle one cell. This is to save the data efficiently, because different circles
% can have different trajectory points
[MRStruct,AdditionalOut.RecoOperators] = op_ReconstructNonCartMRData(MRStruct,[],Settings.NonCartReco);


%% DEBUG: PLOT NonCart Trajectories

if(Settings.Debug.ShowTrajs)
    PlotPartition = ceil(numel(MRStruct.InTraj.GM)/2);
    figure;
    scatter(squeeze(MRStruct.OutTraj.GM(1,1,:)),squeeze(MRStruct.OutTraj.GM(2,1,:)),'b'), hold on 
    CurInd = 1;
    for AngIntNo = 1:MRStruct.Par.nCirc(PlotPartition)
        CurTrajPts = MRStruct.Par.TrajPts(1);
        scatter(squeeze(MRStruct.InTraj.GM{PlotPartition}(1,CurInd+2:CurInd-1+CurTrajPts)), squeeze(MRStruct.InTraj.GM{PlotPartition}(2,CurInd+2:CurInd-1+CurTrajPts)),30,'r')
        plot(squeeze(MRStruct.InTraj.GM{PlotPartition}(1,CurInd:CurInd-1+CurTrajPts)), squeeze(MRStruct.InTraj.GM{PlotPartition}(2,CurInd:CurInd-1+CurTrajPts)),'r')
        scatter(squeeze(MRStruct.InTraj.GM{PlotPartition}(1,CurInd)), squeeze(MRStruct.InTraj.GM{PlotPartition}(2,CurInd)),30,'m')
        scatter(squeeze(MRStruct.InTraj.GM{PlotPartition}(1,CurInd+1)), squeeze(MRStruct.InTraj.GM{PlotPartition}(2,CurInd+1)),	30,'g')
        CurInd = CurInd + CurTrajPts;


    end
    hold off
end


%% Slice Reco

MRStruct = op_SliceReco(MRStruct);



