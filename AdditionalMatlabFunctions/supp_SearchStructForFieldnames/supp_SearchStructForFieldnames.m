function [FoundString, Findings] = supp_SearchStructForFieldnames(Structure,SearchString)
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

FoundString = [];


%% Search Top Layer

FirstLayerFieldnames = fieldnames(Structure); 

FoundInd = 0;
Found_Logi = contains(FirstLayerFieldnames,SearchString);
Occurrences = sum(Found_Logi);
if(Occurrences > 0)
    Findings.FirstLayer = FirstLayerFieldnames(Found_Logi);
    Tmp = cellfun(@(f) Structure.(f), FirstLayerFieldnames(Found_Logi), 'UniformOutput', false);
    Findings.FirstLayer(:,2) = Tmp;
    FoundInd = FoundInd + 1;
    FoundString{FoundInd} = ['FirstLayer: ' num2str(Occurrences)];
end



%% Search First Layer

for ii = 1:numel(FirstLayerFieldnames)
    CurStruct = Structure.(FirstLayerFieldnames{ii});
    if(isstruct(CurStruct))
        SecondLayerFieldnames = fieldnames(Structure.(FirstLayerFieldnames{ii}));
        Found_Logi = contains(SecondLayerFieldnames,SearchString);
        Occurrences = sum(Found_Logi);

        if(Occurrences > 0)
            Findings.(FirstLayerFieldnames{ii}) = SecondLayerFieldnames(Found_Logi);
            Tmp = cellfun(@(f) CurStruct.(f), SecondLayerFieldnames(Found_Logi), 'UniformOutput', false);
            Findings.(FirstLayerFieldnames{ii})(:,2) = Tmp;
            FoundInd = FoundInd + 1;
            FoundString{FoundInd} = ['FirstLayer.' FirstLayerFieldnames{ii} ': ' num2str(Occurrences)];
        end        
        
    end
    
end


