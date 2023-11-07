function [MRStruct,AdditionalOut] = op_CalcNoiseCorrMat(MRStruct,Settings)
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
%                                           CartFFT:        Settings for the non-Cartesian MRSI Reco, see op_ReconstructNonCartMRData.m.
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

if(~isfield_recursive(Settings,'CutAwayFirstNPoints'))
   Settings.CutAwayFirstNPoints = 19;    
end

if(Settings.CutAwayFirstNPoints < 1)
    Settings.CutAwayFirstNPoints = 1;
end





%% Permute, Reshape & Cut Data

% Assume data comes from unlocalized FID sequence
if(numel(size(MRStruct.Data)) > 2)
    MRStruct = op_PermuteMRData(MRStruct,[6 5 7 1 2 3 4]);
    MRStruct = op_ReshapeMRData(MRStruct,[MRStruct.RecoPar.DataSize(1) prod(MRStruct.RecoPar.DataSize(2:end))]);
else
    MRStruct = op_PermuteMRData(MRStruct,[2 1]);
end


% Cut away first N points, they sometimes contain signal
if(Settings.CutAwayFirstNPoints > size(MRStruct.Data,2))
    Settings.CutAwayFirstNPoints = size(MRStruct.Data,2);
end
MRStruct.Data = MRStruct.Data(:,Settings.CutAwayFirstNPoints:end);



%% Calc NoiseCorrMat

if(nargout > 1)
    AdditionalOut.NoiseData = MRStruct;
end

MRStruct.Data = 1/(size(MRStruct.Data,2)) * (MRStruct.Data * MRStruct.Data');
MRStruct.RecoPar.DataSize = size(MRStruct.Data);


%% Postparations
MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);




