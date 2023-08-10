function [MRStruct,AdditionalOut] = io_ReadAverageReshapeImagingGREData(file)
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
% -         file          ...  The Siemens twix-data file of the spiral sequence   
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
MRStruct.DataFile = file;


%% Read Data

% Initial size: [nAngInt*nTempInt x samples x nADCs x nCha x nAvg x nPart x nSlc]
MRStruct.Data = mapVBVD(file);

% ver = MRStruct.Data.image.softwareVersion;
MRStruct.Data = MRStruct.Data.image();
MRStruct = supp_UpdateRecoSteps(MRStruct,struct(),'mapVBVD');

MRStruct.Par = read_ascconv(file);

MRStruct.Data = single(MRStruct.Data);


%% Reshape Data

MRStruct = op_PermuteMRData(MRStruct,[1 3 4 5 7 2 10 11 6 8 9 12 13 14 15]);

% % Assymmetric Echo Zerofilling
% if(MRStruct.RecoPar.DataSize(1) < MRStruct.Par.nFreqEnc*2 )
%     NewDataSize = MRStruct.RecoPar.DataSize; NewDataSize(1) = MRStruct.Par.nFreqEnc*2;
%     MRStruct = op_XFillOrCutData(MRStruct,struct('X',0,'Zerofill_To',NewDataSize,'AppendZerosTo',{{'Beginning'}}));
% end
MRStruct.Data = {MRStruct.Data};

    
end




