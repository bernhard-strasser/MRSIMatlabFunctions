function [MRStruct,MRStructRaw] = op_ReshapeGenericMRData(MRStructRaw,CurDataSet)

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



    %% 0. Preparations

    if(~exist('Settings','var'))
       Settings = struct; 
    end
    MRStruct = struct;
    
    % Find exact name
    CurFields = fieldnames(MRStructRaw.Data);
    FoundField = find(~cellfun(@isempty,regexpi(CurFields,CurDataSet)));
    if(numel(FoundField) > 1)
        fprintf('\nWarning in op_ReshapeGenericMRData: Found several data sets containing name %s',CurDataSet)
    elseif(numel(FoundField) == 0)
        fprintf('\nWarning in op_ReshapeGenericMRData: Found no data set containing name %s',CurDataSet)
        return;
    end
    CurDataSet = CurFields{FoundField(1)};
    
    
    MRStruct.Par = MRStructRaw.Par;
    MRStruct.RecoSteps = MRStructRaw.RecoSteps;
    
    
    %%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape raw Data 1 & Perform averaging                                               %
    % Current Size: {nTotADCs}[nADCPts]                                                    %
    % Target Size:  {nCircles}[nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs x nCha] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Data
    MRStruct.Par.dimnames = {'Lin','Phs','Par','Sli','Col','Cha','Ave','Eco','Rep','Set','Seg','Ida','Idb','Idc','Idd','Ide'};
    MRStruct.Par.DataSize = [...
    MRStructRaw.mdhInfo.(CurDataSet).NLin MRStructRaw.mdhInfo.(CurDataSet).NPhs MRStructRaw.mdhInfo.(CurDataSet).NPar MRStructRaw.mdhInfo.(CurDataSet).NSli ...
    MRStructRaw.mdhInfo.(CurDataSet).NCol MRStructRaw.mdhInfo.(CurDataSet).NCha MRStructRaw.mdhInfo.(CurDataSet).NAve MRStructRaw.mdhInfo.(CurDataSet).NEco ...
    MRStructRaw.mdhInfo.(CurDataSet).NRep MRStructRaw.mdhInfo.(CurDataSet).NSet MRStructRaw.mdhInfo.(CurDataSet).NSeg MRStructRaw.mdhInfo.(CurDataSet).NIda ...
    MRStructRaw.mdhInfo.(CurDataSet).NIdb MRStructRaw.mdhInfo.(CurDataSet).NIdc MRStructRaw.mdhInfo.(CurDataSet).NIdd MRStructRaw.mdhInfo.(CurDataSet).NIde];
    MRStruct.Data = zeros(MRStruct.Par.DataSize,'single'); 

    % Rearrange Data
    for ii = 1:numel(MRStructRaw.Data.(CurDataSet))
        MRStruct.Data(...
        MRStructRaw.mdhInfo.(CurDataSet).Lin(ii),MRStructRaw.mdhInfo.(CurDataSet).Phs(ii),MRStructRaw.mdhInfo.(CurDataSet).Par(ii), ...
        MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,:, ...
        MRStructRaw.mdhInfo.(CurDataSet).Ave(ii),MRStructRaw.mdhInfo.(CurDataSet).Eco(ii),MRStructRaw.mdhInfo.(CurDataSet).Rep(ii), ...
        MRStructRaw.mdhInfo.(CurDataSet).Set(ii),MRStructRaw.mdhInfo.(CurDataSet).Seg(ii),MRStructRaw.mdhInfo.(CurDataSet).Ida(ii), ...
        MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),MRStructRaw.mdhInfo.(CurDataSet).Idc(ii),MRStructRaw.mdhInfo.(CurDataSet).Idd(ii),...
        MRStructRaw.mdhInfo.(CurDataSet).Ide(ii)) ...
        = MRStructRaw.Data.(CurDataSet){ii};
        
        % Automatic averaging of 'set' dimension
%         TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)}(MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),MRStructRaw.mdhInfo.(CurDataSet).Seg(ii),MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Ida(ii),:) = ...
%         TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)}(MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),MRStructRaw.mdhInfo.(CurDataSet).Seg(ii),MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Ida(ii),:) + CurData; 
    end
%     TmpData2 = cellfun(@(x) x/MRStructRaw.mdhInfo.(CurDataSet).NSet,TmpData2,'UniformOutput',false);     % For averaging
%     clear TmpData1; 


    MRStruct.Par.DataSize = size(MRStruct.Data);     % Generalize later for more than 1 partition!
    MRStruct.Par.dimnames = MRStruct.Par.dimnames(1:numel(MRStruct.Par.DataSize));
    MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,CurDataSet);
    MRStructRaw.Data = rmfield(MRStructRaw.Data,CurDataSet);
    
    
    