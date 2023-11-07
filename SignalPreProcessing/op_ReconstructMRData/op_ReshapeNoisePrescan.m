function [MRStruct,MRStructRaw] = op_ReshapeNoisePrescan(MRStructRaw)    

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
    % [MRStruct, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
    %
    % Input: 
    % -         DataFile          ...  The Siemens twix-data file of the spiral sequence   
    % -         TrajectoryFile    ...  The trajectory file containing information about the k-space spiral trajectory.
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape raw Data 1 & Perform averaging                                               %
    % Current Size: {nTotADCs}[nADCPts]                                                    %
    % Target Size:  {nCircles}[nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs x nCha] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Preparations
    Settings = struct;
    CurDataSet = 'NOISEADJSCAN';

    if(~isfield(MRStructRaw.Data,CurDataSet))
        fprintf('\nWarning in op_ReshapeNoisePrescan: Found no data set containing name %s',CurDataSet)        
        MRStruct = struct;
        return;
    end
    
    %%
    
    
    % Rearrange Data
    MRStruct.Data = cat(1,MRStructRaw.Data.NOISEADJSCAN{:});
    MRStruct.Par = MRStructRaw.Par;
    
    Fieldies = fieldnames(MRStruct.Par);
    for SearchFor = {'VoI','FoV','Pos','NormalVector','TEs','Slice','InPlaneRotation'}
        MRStruct.Par = supp_SetAllFields(MRStruct.Par,Fieldies(~cellfun(@isempty,regexpi(Fieldies,SearchFor{1}))),NaN);
    end
    if(isfield_recursive(MRStruct.Par,'WipMemBlockInterpretation.Prescan.NOISEADJSCAN.Dwelltime'))
        MRStruct.Par.Dwelltimes = MRStruct.Par.WipMemBlockInterpretation.Prescan.NOISEADJSCAN.Dwelltime;
    else    % Just guess
        MRStruct.Par.Dwelltimes = 1E4;
    end
    for SearchFor = {'Enc','vecSize'}
        MRStruct.Par = supp_SetAllFields(MRStruct.Par,Fieldies(~cellfun(@isempty,regexpi(Fieldies,SearchFor{1}))),1);
    end
    MRStruct.Par.nPtsPerADC = MRStructRaw.mdhInfo.NOISEADJSCAN.Col;
    MRStruct.Par.nADCs = MRStructRaw.mdhInfo.NOISEADJSCAN.NLin;
    
    
    MRStruct.Par.DataSize = size(MRStruct.Data);     % Generalize later for more than 1 partition!
    MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,CurDataSet);
    MRStructRaw.Data = rmfield(MRStructRaw.Data,CurDataSet);
    
    
%%

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);

    
