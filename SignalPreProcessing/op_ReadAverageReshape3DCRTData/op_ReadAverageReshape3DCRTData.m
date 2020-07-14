function [MRStruct,MRStruct_refscan] = op_ReadAverageReshape3DCRTData(MRStruct, Settings)
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

    if(~isfield(Settings,'ProduceNoiseThroughAvgs_flag'))
        Settings.ProduceNoiseThroughAvgs_flag = true;
    end
    if(~isfield(Settings,'ReadInAdditionalDataSets'))
        Settings.ReadInAdditionalDataSets = {};
    end

    % Time total
    TicTotal = tic;
    fprintf('\n\nRead and Reshape data\t\t\t...')

    %% Read DataStructure

    fprintf('\n\tReading data structure\t\t\t...\n')
    ticcyRead = tic;
    Tmp = mapVBVD(MRStruct.DataFile);
    MRStruct.MapVBVDObj = Tmp.image;

    MRStruct_refscan = MRStruct;
    MRStruct_refscan.MapVBVDObj = Tmp.refscan;
    fprintf('\n\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyRead))


    %% Read Other Datasets

    % for CurSet = 1:numel(Settings.ReadInAdditionalDataSets)
    %     ticcyReadOthers = tic;
    %     fprintf('\n\tReading data %s\t\t\t...\n',Settings.ReadInAdditionalDataSets{CurSet})
    %     MRStruct.(Settings.ReadInAdditionalDataSets{CurSet}) = Tmp.(Settings.ReadInAdditionalDataSets{CurSet})();
    %     fprintf('\n\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyReadOthers))
    % end
    MRStruct.NoiseData = Tmp.noise();
    MRStruct_refscan.NoiseData = Tmp.noise();
    clear Tmp;


    %% Read Parameters

    fprintf('\n\tReading Parameters\t\t\t...\n')
    ticcyReadPars = tic;
    MRStruct = io_Read3DCRTPars(MRStruct,struct());
    MRStruct_refscan = io_Read3DCRTPars(MRStruct_refscan,struct());
    fprintf('\n\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyReadPars))


    %% Read & Reshape Data 
    ticcyReadAverageReshape = tic;
    fprintf('\n\tReading, Averaging & Reshaping data\t...\n')
    MRStruct = ReadAndReshapeCRTData(MRStruct);
    MRStruct_refscan = ReadAndReshapeCRTData(MRStruct_refscan);    
    fprintf('\n\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyReadAverageReshape))


    %% Postparations

    % Output.Data = MRStruct.Data;
    % if(isfield(MRStruct,'NoiseData') && ~isempty(MRStruct.NoiseData))
    %     Output.NoiseData = MRStruct.NoiseData;
    % end
    % clear MRStruct
    % Output.MRStruct.Par = MRStruct.Par;

    %% The End

    MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);
    MRStruct_refscan = supp_UpdateRecoSteps(MRStruct_refscan,Settings);
    
    
    fprintf('\nTotal Read, Average, Reshape\t\t...\ttook\t%10.6f seconds',toc(TicTotal))
end


%%
function MRStruct = ReadAndReshapeCRTData(MRStruct)

    MRStruct.MapVBVDObj.flagAverageDim = false([1 16]); MRStruct.MapVBVDObj.flagAverageDim(10) = true;          % Does that save memory?

    % Initialize Data
    TmpData2 = cell([1 MRStruct.Par.nAngInts]); 

    for CurCrcl = 1:MRStruct.Par.nAngInts
        TmpData2{CurCrcl} = single(zeros([MRStruct.Par.total_channel_no_reco MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC(CurCrcl) MRStruct.Par.nADCsPerAngInt(CurCrcl) ])); 
    end

    % Check if enough memory is available to read everything in at once
    [~,CurMemFree] = memused_linux(1);
    MemNecessary = sum(cellfun(@numel,TmpData2))*MRStruct.Par.nAverages*2*4/2^20;
    if(1.3*MemNecessary < CurMemFree)

        TmpData1 = single(MRStruct.MapVBVDObj.unsorted());
        MRStruct = supp_UpdateRecoSteps(MRStruct,struct(),'mapVBVD');


        % Reshape raw Data 1 & Perform averaging

        % Get data. Perform memory-efficient averaging
        % TODO: Create NoiseData From Averages
        % Current Size {nCircles}[nCha x nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs]
        for ii = 1:size(TmpData1,3)
            CurData = reshape(transpose(TmpData1(1:MRStruct.MapVBVDObj.Col(ii),:,ii)),size( TmpData2{MRStruct.MapVBVDObj.Lin(ii)}(:,1,1,1,:,1)));
            TmpData2{MRStruct.MapVBVDObj.Lin(ii)}(:,MRStruct.MapVBVDObj.Idb(ii),MRStruct.MapVBVDObj.Seg(ii),MRStruct.MapVBVDObj.Sli(ii),:,MRStruct.MapVBVDObj.Ida(ii)) = TmpData2{MRStruct.MapVBVDObj.Lin(ii)}(:,MRStruct.MapVBVDObj.Idb(ii),MRStruct.MapVBVDObj.Seg(ii),MRStruct.MapVBVDObj.Sli(ii),:,MRStruct.MapVBVDObj.Ida(ii)) + CurData; 
        end
        TmpData2 = cellfun(@(x) x/max(MRStruct.MapVBVDObj.Set),TmpData2,'UniformOutput',false);     % For averaging
        clear TmpData1; 
        
    else
        for CurCrcl = 1:MRStruct.Par.nAngInts
            TmpData2{CurCrcl} = sum(MRStruct.MapVBVDObj(1:MRStruct.Par.nPtsPerADC(CurCrcl),:,CurCrcl,:,:,:,:,:,:,:,:,1:MRStruct.Par.nADCsPerAngInt(CurCrcl), 1:MRStruct.Par.nTempIntsPerAngInt(CurCrcl),MRStruct.MapVBVDObj.NIdc,MRStruct.Par.nTempIntsPerAngInt(CurCrcl)+1),10);   % MRStruct.MapVBVDObj.NIdc ok? Was 181
            TmpData2{CurCrcl} = permute(TmpData2{CurCrcl},[2 13 4 5 1 12    3 6 7 8 9 10 11 14 15 16]);
        end
    end




    % Reshape Data to Output Size {MRStruct.Par.nAngInts}[nTrajPoints x nPart x nSlc x vecSize x nCha]
    MRStruct.Data = []; MRStruct = rmfield(MRStruct,'Data');
    for Curkz = 1:MRStruct.Par.nPartEnc
        MRStruct.Data{Curkz} = single(zeros([MRStruct.Par.TrajPts MRStruct.Par.nAngInts 1 MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]));
    end
    % for CurCrcl = 1:MRStruct.Par.nAngInts
    %     MRStruct.Data{CurCrcl} = zeros([MRStruct.Par.TrajPts MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]);
    % end
    for CurCrcl = 1:MRStruct.Par.nAngInts
        TmpData3 = reshape(TmpData2{CurCrcl},[MRStruct.Par.total_channel_no_reco MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC(CurCrcl)*MRStruct.Par.nADCsPerAngInt(CurCrcl) ]); 
        TmpData3 = TmpData3(:,:,:,:,1:MRStruct.Par.UsefulADCPtsPerAngInt(CurCrcl));
        TmpData3 = reshape(TmpData3,[MRStruct.Par.total_channel_no_reco MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.TrajPts MRStruct.Par.vecSize/MRStruct.Par.nTempIntsPerAngInt(CurCrcl)]); 
        TmpData3 = permute(TmpData3,[5 3 4 2 6 1]); % New size : [nTrajPoints x nPart x nSlc x vecSize x nTempIntsPerAngInt x nCha]
        TmpData3 = reshape(TmpData3,[MRStruct.Par.TrajPts 1 MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]);
    %     MRStruct.Data{CurCrcl} = TmpData3;
        for Curkz = 1:MRStruct.Par.nPartEnc
            MRStruct.Data{Curkz}(:,CurCrcl,:,:,:,:) = TmpData3(:,:,Curkz,:,:,:);
        end
    end


    MRStruct.Par.DataSize = size(MRStruct.Data{1});     % Generalize later for more than 1 partition!
    MRStruct = rmfield(MRStruct,'MapVBVDObj');



end
