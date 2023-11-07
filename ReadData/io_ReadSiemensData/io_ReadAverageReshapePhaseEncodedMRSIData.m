function [MRStruct,RefStruct,NoiseStruct] = io_ReadAverageReshapePhaseEncodedMRSIData(file)
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


%% Read Data

MRStructRaw = io_ReadSiemensData(file);


%% Create WeightedAveragingMap

for CurSet2 = transpose(fieldnames(MRStructRaw.mdhInfo))
    CurSet = CurSet2{:};

    WeightedAveragingMap.(CurSet) = zeros([MRStructRaw.mdhInfo.(CurSet).NLin MRStructRaw.mdhInfo.(CurSet).NSeg MRStructRaw.mdhInfo.(CurSet).NPar]); 
    for ii = 1:numel(MRStructRaw.mdhInfo.(CurSet).Ave)
        WeightedAveragingMap.(CurSet)(MRStructRaw.mdhInfo.(CurSet).Lin(ii),MRStructRaw.mdhInfo.(CurSet).Seg(ii),MRStructRaw.mdhInfo.(CurSet).Par(ii)) = MRStructRaw.mdhInfo.(CurSet).Ave(ii); 
    end
end



%% Reshape Data


[MRStruct,MRStructRaw] = op_ReshapeGenericMRData(MRStructRaw,'ONLINE');
MRStruct.WeightedAveragingMap = WeightedAveragingMap.ONLINE;
[RefStruct,MRStructRaw] = op_ReshapeGenericMRData(MRStructRaw,'PATREF*');
if(numel(fieldnames(RefStruct)) > 0)
    MRStruct.WeightedAveragingMap = WeightedAveragingMap.PATREFSCAN;
end
NoiseStruct = op_ReshapeNoisePrescan(MRStructRaw);




%% Correct RefStruct fields

if(numel(fieldnames(RefStruct)) > 0)
    RefStruct.Par.AssumedSequence = 'Imaging_GRE';
    RefStruct.Par.nFreqEnc = size(RefStruct.Data,5)/2;
    RefStruct.Par.nPhasEnc = size(RefStruct.Data,1);
    RefStruct.Par.nFreqEnc_FinalMatrix = size(RefStruct.Data,5);
    RefStruct.Par.nPhasEnc_FinalMatrix = size(RefStruct.Data,1);
    RefStruct.Par.SpatialSpectralEncoding_flag = false;
end
MRStruct.Par.SpatialSpectralEncoding_flag = false;


%% Reshape MRI Data


if(regexp(MRStruct.Par.ScannerVersion,'syngo MR B'))
    RefStruct = op_PermuteMRData(RefStruct,[5 1 3 4 2 6 10 7 8 9]);
%     RefStruct = op_PermuteMRData(RefStruct,[1 3 4 5 7 2 10 6 8 9 11 12 13 14 15 16]);     % For MapVBVD
%     MRStruct = op_PermuteMRData(MRStruct,[3 7 4 5 1 2 6 8 9 10 11 12 13 14 15 16]);
else
    MRStruct = op_PermuteMRData(MRStruct,[11 1 3 4 5 6 7 2 8 9 10]);
end


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);
RefStruct = supp_UpdateRecoSteps(RefStruct,Settings);
NoiseStruct = supp_UpdateRecoSteps(NoiseStruct,Settings);



