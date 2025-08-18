function [FittingResultMaps,FittingResultSpectra,AdditionalOut] = op_PerformLCModelFitting(MRStruct,FittingMask,Paths,LCModelControlInfo,Settings)
%
% op_ReadAndRecoBorjanSpiralData Read and reconstruct data from Borjan Gagoski's Spiral MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in Spiral MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         SpiralDataFile          ...  The Siemens twix-data file of the spiral sequence   
% -         SpiralTrajectoryFile    ...  The trajectory file containing information about the k-space spiral trajectory.
% -         Settings                ...  Struct with fields
%                                           Debug: The Debug-settings. Subfields:
%                                                           ShowTrajs:  If you want to plot the spiral and Cartesian trajectories to get a feeling how you measured...
%                                           io_ReadSpiralPars:  Settings for the function "io_ReadSpiralPars", see io_ReadSpiralPars.m
%                                           ReadInTraj:         Settings for reading the trajectory, see io_ReadSpiralTraj.m.
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
%                                           InTraj:     The (spiral) trajectory with which the data were measured.
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
if(~isfield(Settings,'NoOfCPUCores'))
    Settings.NoOfCPUCores = min(16,sum(FittingMask(:)));
end

if(~isfield(Settings,'AllowNonEmptyDirs'))
    Settings.AllowNonEmptyDirs = false;
end

if(~exist('LCModelControlInfo','var'))
    LCModelControlInfo = 0;
end

AdditionalOut = struct;


%% Create & Check Paths

if(~exist(Paths.LCM_Path,'file'))
    error('From op_PerformLCModelFitting: Please provide correct path to LCModel') 
end

Tmp = dir(Paths.OutputDir);
if(numel(Tmp) > 0)  % Folder exists 
    if(numel(Tmp) > 2) % Folder exists and has more content than just  '.' and '..' (i.e. it is not empty)
        if(~Settings.AllowNonEmptyDirs)
            error('From op_PerformLCModelFitting: Output directory is not empty.')
        else
            warning('From op_PerformLCModelFitting: Output directory is not empty.')
        end
    end
else
    mkdir(Paths.OutputDir);
end

Tmp = dir(Paths.TmpDir);
if(numel(Tmp) > 0)  % Folder exists 
    if(numel(Tmp) > 2) % Folder exists and has more content than just  '.' and '..' (i.e. it is not empty)
        if(~Settings.AllowNonEmptyDirs)
            error('From op_PerformLCModelFitting: TmpDir directory is not empty.')
        else
            warning('From op_PerformLCModelFitting: TmpDir directory is not empty.')
            delete([Paths.TmpDir '/lcm_process_core_*.sh'])
            delete([Paths.TmpDir '/RunLCModel.sh'])
        end
    end
else
    mkdir(Paths.TmpDir);
end


%% Write LCModel Files

	
	
	PathsForWriteLCMFiles.out_dir = Paths.OutputDir;
	PathsForWriteLCMFiles.basis_file = {Paths.BasisSetFile};                            % Path to Basis set that should be used ('.basis')
	PathsForWriteLCMFiles.LCM_ProgramPath = Paths.LCM_Path;                         % Path of the LCModel program
	PathsForWriteLCMFiles.batchdir = Paths.TmpDir;
	MetaInfo.DimNames = {'x','y','z'};
    
    if(isfield(MRStruct.RecoPar,'PatientName'))
        MetaInfo.pat_name = MRStruct.RecoPar.PatientName;                                % Name of Patient which determines the naming of the output spectra
    else
        MetaInfo.pat_name = 'UnknownName';
    end
	MetaInfo.LarmorFreq = MRStruct.RecoPar.LarmorFreq;
	MetaInfo.dwelltime = MRStruct.RecoPar.Dwelltimes(1);
    
	if(isfield(Paths,'priors_dir'))
        fprintf('\nUse priors from OFF spectra.\nWARNING: THIS OPTION WORKS ONLY FOR USING MEGA-OFF PRIOR KNOWLEDGE FOR MEGA-DIFF-FITTING\n' )
		PathsForWriteLCMFiles.priors_dir = Paths.priors_dir;
	end	
    
	Data.csi = MRStruct.Data;
	Write_LCM_files(Data,PathsForWriteLCMFiles,MetaInfo,LCModelControlInfo,FittingMask,Settings.NoOfCPUCores,true) 




%% Run LCModel

ticcy = tic;
fprintf('\n\nStarting LCModel fitting')
system(['bash ' Paths.TmpDir '/RunLCModel.sh']);
fprintf('\n\nFinished LCModel fitting within %f s',toc(ticcy))

% Remove unnecessary files

delete([Paths.OutputDir '/' MetaInfo.pat_name '_*.RAW'])
delete([Paths.OutputDir '/' MetaInfo.pat_name '_*.control'])
delete([Paths.TmpDir '/lcm_process_core_*.sh'])
delete([Paths.TmpDir '/Core_*_Finished.txt'])
delete([Paths.TmpDir '/RunLCModel.sh'])
TmpDirr = dir(Paths.TmpDir);
if(numel(TmpDirr) == 2)
    rmdir(Paths.TmpDir)
end
clear TmpDirr 

%% Read Data


FittingResultMaps = io_ReadAllTableFiles(Paths.OutputDir,MRStruct);
FittingResultSpectra = read_CoordFolder([Paths.OutputDir '/CoordFiles'],MRStruct,false);




