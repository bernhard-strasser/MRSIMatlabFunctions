function [MRStruct,AdditionalOut] = op_ReconstructCartMRSI(MRStruct,Settings)
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
if(~isfield_recursive(Settings,'Averaging.UndoWeightedAveraging_flag'))
    Settings.Averaging.UndoWeightedAveraging_flag = true;
end


%% Sum Averages

[MRStruct] = op_AverageMRData(MRStruct,Settings.Averaging);



%% Fourier Transform Data

if(MRStruct.RecoPar.nPartEnc > 1)
    Settings.zCartFFT.ConjFlag = false;    
    Settings.zCartFFT.Ifft_flag = false;    
    Settings.zCartFFT.FlipDim_flag = false;
    Settings.zCartFFT.ApplyAlongDims = [3];
    MRStruct = op_FFTOfMRIData_v2(MRStruct,Settings.zCartFFT);        
end

Settings.CartFFT.ConjFlag = true;    
Settings.CartFFT.Ifft_flag = false;    
Settings.CartFFT.FlipDim_flag = true;
Settings.CartFFT.FlipDim = 1;
Settings.CartFFT.ApplyAlongDims = [1 2];
MRStruct = op_FFTOfMRIData_v2(MRStruct,Settings.CartFFT);
clear Settings;

% MRStruct = op_FFTOfMRIData_v2(MRStruct,struct('ApplyAlongDims',[1 2 3],'ConjFlag',1,'FlipDim',1));


%% Slice Reco

MRStruct = op_SliceReco(MRStruct);




