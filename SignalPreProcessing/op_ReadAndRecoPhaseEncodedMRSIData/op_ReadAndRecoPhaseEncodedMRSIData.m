function [MRStruct,RefScan, NoiseScan] = op_ReadAndRecoPhaseEncodedMRSIData(DataFile,Settings,AdditionalInput)
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
if(~exist('AdditionalInput','var'))
    AdditionalInput = struct;
end

if(~isfield_recursive(Settings,'CartFFT.ConjFlag'))
   Settings.CartFFT.ConjFlag = true;    
end
if(~isfield_recursive(Settings,'CartFFT.Ifft_flag'))
   Settings.CartFFT.Ifft_flag = false;    
end
if(~isfield_recursive(Settings,'CartFFT.FlipDim_flag'))
    Settings.CartFFT.FlipDim_flag = true;
end
if(~isfield_recursive(Settings,'CartFFT.FlipDim'))
    Settings.CartFFT.FlipDim = 1;
end
if(~isfield_recursive(Settings,'CartFFT.ApplyAlongDims'))
    Settings.CartFFT.ApplyAlongDims = [1 2];
end
if(~isfield_recursive(Settings,'zCartFFT.ConjFlag'))
   Settings.zCartFFT.ConjFlag = false;    
end
if(~isfield_recursive(Settings,'zCartFFT.Ifft_flag'))
   Settings.zCartFFT.Ifft_flag = false;    
end
if(~isfield_recursive(Settings,'zCartFFT.FlipDim_flag'))
    Settings.zCartFFT.FlipDim_flag = false;
end
if(~isfield_recursive(Settings,'zCartFFT.ApplyAlongDims'))
    Settings.zCartFFT.ApplyAlongDims = [3];
end

if(~exist('Settings','var'))
    Settings = struct;
end


% Create MRStruct
if(isstruct(DataFile) && isfield(DataFile,'DataFile'))
    MRStruct = DataFile;
else
    MRStruct.DataFile = DataFile;
end






%% Read, Average & Reshape Data

[MRStruct, RefScan, NoiseScan] = io_ReadAndReshapeSiemensData(MRStruct.DataFile);

% MRStruct = RefScan;



%% Noise-decorrelate data

if(isfield(NoiseScan,'Data') && ~isfield(AdditionalInput,'NoiseCorrMatStruct'))
    % Create AdditionalInput.Noise
end
    
if(isfield(AdditionalInput,'NoiseCorrMatStruct'))
   
    MRStruct = op_PerformNoiseDecorrelation(MRStruct,AdditionalInput.NoiseCorrMatStruct);
    if(isfield(RefScan,'Data'))
        RefScan = op_PerformNoiseDecorrelation(RefScan,AdditionalInput.NoiseCorrMatStruct);
    end
    
end


%% Reconstruct Cart Data

if(~MRStruct.Par.dicom_flag)
    
    % Do z-fft before conj that is sometimes done in in-plane FFT
    if(MRStruct.Par.nPartEnc > 1)
        MRStruct = op_FFTOfMRIData_v2(MRStruct,Settings.zCartFFT);        
    end
    
    MRStruct = op_FFTOfMRIData_v2(MRStruct,Settings.CartFFT);
end


%% Slice Reco

MRStruct = op_SliceReco(MRStruct);



