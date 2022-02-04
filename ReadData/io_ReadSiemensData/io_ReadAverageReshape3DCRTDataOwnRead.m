function [MRStruct,MRStruct_refscan,MRStruct_Noise] = io_ReadAverageReshape3DCRTDataOwnRead(file, Settings)
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
    MRStruct_Noise = [];
    % Time total
    TicTotal = tic;
    fprintf('\n\nRead and Reshape data\t\t\t\t...')

    %% Read Data

    MRStructRaw = io_ReadSiemensData(file);
    
    
    

    %% Read Parameters

    fprintf('\tReading Parameters\t\t\t...')
    ticcyReadPars = tic;
    MRStruct = io_Read3DCRTParsOwnRead(MRStructRaw,'ONLINE'); 
    MRStruct_refscan = io_Read3DCRTParsOwnRead(MRStructRaw,'PATREFSCAN');
    MRStruct.Par.dicom_flag = false;
    MRStruct_refscan.Par.dicom_flag = false;
    fprintf('\n\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyReadPars))


    %% Reshape Data 
    ticcyReadAverageReshape = tic;
    fprintf('\n\tReading, Averaging & Reshaping data\t...\n')
    [MRStruct,MRStructRaw] = Reshape3DCRTDataOwnRead(MRStructRaw,MRStruct,'ONLINE');
    if(nargout > 1)
    	[MRStruct_refscan,MRStructRaw] = Reshape3DCRTDataOwnRead(MRStructRaw,MRStruct_refscan,'PATREFSCAN');
    end
    if(nargout > 2)
    	MRStruct_Noise = ReshapeNoisePrescan(MRStructRaw);
    end
    clear MRStructRaw;
    fprintf('\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyReadAverageReshape))


    %% Postparations

    % Output.Data = MRStruct.Data;
    % if(isfield(MRStruct,'NoiseData') && ~isempty(MRStruct.NoiseData))
    %     Output.NoiseData = MRStruct.NoiseData;
    % end
    % clear MRStruct
    % Output.MRStruct.Par = MRStruct.Par;

    %% The End

    MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);
    if(nargout > 1)    
    	MRStruct_refscan = supp_UpdateRecoSteps(MRStruct_refscan,Settings);
    end
    
    fprintf('\nTotal Read, Average, Reshape\t\t\t...\ttook\t%10.6f seconds',toc(TicTotal))
end


%%
function [MRStruct,MRStructRaw] = Reshape3DCRTDataOwnRead(MRStructRaw,MRStruct,CurDataSet)

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape raw Data 1 & Perform averaging                                               %
    % Current Size: {nTotADCs}[nADCPts]                                                    %
    % Target Size:  {nCircles}[nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs x nCha] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Data
    TmpData2 = cell([1 MRStruct.Par.nAngInts]); 
    for CurCrcl = 1:MRStruct.Par.nAngInts
        TmpData2{CurCrcl} = single(zeros([MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC(CurCrcl) MRStruct.Par.nADCsPerAngInt(CurCrcl) MRStruct.Par.total_channel_no_reco])); 
    end
    % Rearrange Data
    for ii = 1:numel(MRStructRaw.Data.(CurDataSet))
            
        CurCrcl = MRStructRaw.mdhInfo.(CurDataSet).Lin(ii);
        
        % Remark: The 3rd circle might have partitions (2,3,...,14) saved in .Seg(ii). However, we don't want to zero-fill the other data here,
        % so our TmpData2 only has 13 partitions. Therefore using .Seg(ii) will result in an error, and we need to use the numbers to 
        % (1,2,...,13), which is saved in CurSegToSave
        CurSegToSave = MRStructRaw.mdhInfo.(CurDataSet).Seg(ii) - (MRStruct.Par.nPartEnc-size(TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)},2))/2;
        CurData = reshape(MRStructRaw.Data.(CurDataSet){ii},[1 1 1 MRStruct.Par.nPtsPerADC( CurCrcl) 1 MRStruct.Par.total_channel_no_reco]);
        TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)}(MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),CurSegToSave,MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Ida(ii),:) = ...
        TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)}(MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),CurSegToSave,MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Ida(ii),:) + CurData; 
    end
    TmpData2 = cellfun(@(x) x/MRStructRaw.mdhInfo.(CurDataSet).NSet,TmpData2,'UniformOutput',false);     % For averaging
    clear TmpData1; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape raw Data 2                                                                     %
    % Current Size: {nCircles}[nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs x nCha]   %
    % Target Size:  {nCircles}[nTrajPoints x 1 x nPart x nSlc x vecSize x nCha] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MRStruct.Data = []; MRStruct = rmfield(MRStruct,'Data');
    for CurCrcl = 1:MRStruct.Par.nAngInts
        MRStruct.Data{CurCrcl} = zeros([MRStruct.Par.TrajPts(CurCrcl) 1 MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco],'single');
    end
    for CurCrcl = 1:MRStruct.Par.nAngInts
        TmpData3 = reshape(TmpData2{CurCrcl},[MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC(CurCrcl)*MRStruct.Par.nADCsPerAngInt(CurCrcl) MRStruct.Par.total_channel_no_reco]); 
        TmpData3 = TmpData3(:,:,:,1:MRStruct.Par.UsefulADCPtsPerAngInt(CurCrcl),:);
        TmpData3 = reshape(TmpData3,[MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.TrajPts(CurCrcl) MRStruct.Par.vecSize/MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.total_channel_no_reco] ); 
        TmpData3 = permute(TmpData3,[4 2 3 1 5 6]); % Permute from [nCha x nTempIntsPerAngInt x nPart x nSlc x nTrajPoints x vecSize/nTempIntsPerAngInt] to [nTrajPoints x nPart x nSlc x vecSize/nTempIntsPerAngInt x nTempIntsPerAngInt x nCha]
        TmpData3 = reshape(TmpData3,[MRStruct.Par.TrajPts(CurCrcl) 1 MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]); % Merge vecSize/nTempIntsPerAngInt x nTempIntsPerAngInt to vecSize
        MRStruct.Data{CurCrcl} = TmpData3;
%         for Curkz = 1:MRStruct.Par.nPartEnc
%             MRStruct.Data{Curkz}(:,CurCrcl,:,:,:,:) = TmpData3(:,:,Curkz,:,:,:);
%         end
    end


    MRStruct.Par.DataSize = cellfun(@size,MRStruct.Data,'uni',false);     % Generalize later for more than 1 partition!
    MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,CurDataSet);
    MRStructRaw.Data = rmfield(MRStructRaw.Data,CurDataSet);
    
    
    
end


function [MRStruct,MRStructRaw] = ReshapeNoisePrescan(MRStructRaw)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape raw Data 1 & Perform averaging                                               %
    % Current Size: {nTotADCs}[nADCPts]                                                    %
    % Target Size:  {nCircles}[nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs x nCha] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CurDataSet = 'NOISEADJSCAN';
    
    % Rearrange Data
    MRStruct.Data = cat(1,MRStructRaw.Data.NOISEADJSCAN{:});
    MRStruct.Par = MRStructRaw.Par;
    
    Fieldies = fieldnames(MRStruct.Par);
    for SearchFor = {'Dwelltimes','VoI','FoV','Pos','NormalVector','TEs','Slice','InPlaneRotation'}
        MRStruct.Par = supp_SetAllFields(MRStruct.Par,Fieldies(~cellfun(@isempty,regexpi(Fieldies,SearchFor{1}))),NaN);
    end
    for SearchFor = {'Enc','vecSize'}
        MRStruct.Par = supp_SetAllFields(MRStruct.Par,Fieldies(~cellfun(@isempty,regexpi(Fieldies,SearchFor{1}))),1);
    end
    MRStruct.Par.nPtsPerADC = MRStructRaw.mdhInfo.NOISEADJSCAN.Col;
    MRStruct.Par.nADCs = MRStructRaw.mdhInfo.NOISEADJSCAN.NLin;
    
    
    MRStruct.Par.DataSize = size(MRStruct.Data);     % Generalize later for more than 1 partition!
    MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,CurDataSet);
    MRStructRaw.Data = rmfield(MRStructRaw.Data,CurDataSet);
    
    
    
end




function [MRStruct] = io_Read3DCRTParsOwnRead(MRStructRaw,CurDataset)

    % Copy all fields
    for CurField = transpose(fieldnames(MRStructRaw))
        if(strcmpi(CurField{1},'Data') || strcmpi(CurField{1},'mdhInfo'))
            continue;
        end
        MRStruct.(CurField{1}) = MRStructRaw.(CurField{1});
    end

    MRStruct.Par.nAngInts = MRStructRaw.mdhInfo.PATREFSCAN.NLin;
    MRStruct.Par.ReadoutOSFactor = MRStructRaw.Hdr.Dicom.flReadoutOSFactor;
    MRStruct.Par.Dwelltimes = MRStruct.Par.Dwelltimes*MRStruct.Par.ReadoutOSFactor/2;
    
    MRStruct.Par.nTempIntsPerAngInt = zeros([1 MRStruct.Par.nAngInts]); 
    MRStruct.Par.nPartEncsPerAngInt = zeros([1 MRStruct.Par.nAngInts]); 
    MRStruct.Par.nADCsPerAngInt = zeros([1 MRStruct.Par.nAngInts]);
    MRStruct.Par.nPtsPerADC = zeros([1 MRStruct.Par.nAngInts]);
    MRStruct.Par.TrajPts = zeros([1 MRStruct.Par.nAngInts]);
    for CurCrcl = 1:MRStruct.Par.nAngInts
        MRStruct.Par.nTempIntsPerAngInt(CurCrcl) = max(MRStructRaw.mdhInfo.(CurDataset).Idb(MRStructRaw.mdhInfo.(CurDataset).Lin == CurCrcl));        % Only works for 2D
        MRStruct.Par.nPartEncsPerAngInt(CurCrcl) = numel(unique(MRStructRaw.mdhInfo.(CurDataset).Seg(MRStructRaw.mdhInfo.(CurDataset).Lin == CurCrcl)));    % e.g. from (-3,-2,-1,0,1,2,3) --> (0,1,2,3,4,5,6,7). However, mapVBVD adds +1                                                                                            
                                                                                                % so it should be (1,2,3,4,5,6,7,8).
        MRStruct.Par.nADCsPerAngInt(CurCrcl) = max(MRStructRaw.mdhInfo.(CurDataset).Ida(MRStructRaw.mdhInfo.(CurDataset).Lin == CurCrcl));     
        MRStruct.Par.nPtsPerADC(CurCrcl) = max(MRStructRaw.mdhInfo.(CurDataset).Col(MRStructRaw.mdhInfo.(CurDataset).Lin == CurCrcl));
        MRStruct.Par.TrajPts(CurCrcl) = MRStruct.Par.ReadoutOSFactor*max(MRStructRaw.mdhInfo.(CurDataset).iceParam(6,MRStructRaw.mdhInfo.(CurDataset).Lin == CurCrcl));

    end
        
    if(any(MRStruct.Par.TrajPts == 0))  % For Lukas's sequence, he writes the TrajPts into Idc 
        MRStruct.Par.TrajPts = repmat(MRStruct.Par.ReadoutOSFactor*(MRStructRaw.mdhInfo.(CurDataset).Idc(1)-1),[1 MRStruct.Par.nAngInts]); 
    end

    % Define which angular interleaves are measured for each kz (matrix of nPartEncmax x nCircles)
%     MRStruct.Par.nPartEncsPerAngInt = [9 9 9 7 7 7 7 5 5 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1]; % For simulating 3D
    dummy1 = max(MRStruct.Par.nPartEncsPerAngInt);
    dummy2 = ceil(dummy1/2);
    for kz = 1:max(MRStruct.Par.nPartEncsPerAngInt)
        MRStruct.Par.AngIntsPerPartEnc(kz,:) = MRStruct.Par.nPartEncsPerAngInt >= dummy1-(     (kz - (kz>dummy2)*2*mod(kz,dummy2))        -1)*2;
    end

    % Hack -- Delete me later
    if(MRStruct.Par.TrajPts == 0)
        MRStruct.Par.TrajPts = [4 10 18 24 30 36 45 48 60 60 72 80 80 90 120 120 120 120 120 144 144 144 144 180 180 180 180 180 180 216 216 216] * 2;
    end



    MRStruct.Par.RewPts = 0;
    MRStruct.Par.TrajTotPts = MRStruct.Par.TrajPts + MRStruct.Par.RewPts;



    MRStruct.Par.VTI_Flag = ~all(MRStruct.Par.nTempIntsPerAngInt == MRStruct.Par.nTempIntsPerAngInt(1));



    % For getting vecSize, it depends on whether we read in ONLINE or refscan data
    if(strcmpi(CurDataset,'ONLINE'))
         MRStruct.Par.vecSize = MRStruct.Hdr.Dicom.alICEProgramPara(7);  
    else
        MRStruct.Par.vecSize = MRStruct.Hdr.Dicom.alICEProgramPara(8);
        if(any(mod(MRStruct.Par.vecSize,MRStruct.Par.nTempIntsPerAngInt)))
            MRStruct.Par.vecSize = MRStruct.Hdr.Dicom.alICEProgramPara(8) * max(MRStruct.Par.nTempIntsPerAngInt);    % In one version of my sequence we need to do that...
        end

    end
    if(MRStruct.Par.vecSize <= 0)

        MRStruct.Par.vecSize = round(MRStruct.Par.nPtsPerADC(1)*MRStruct.Par.nADCsPerAngInt(1)/MRStruct.Par.TrajPts(1)*MRStruct.Par.nTempIntsPerAngInt(1)-0.5);

        % vecSize must be dividable by least common multiple of all TIs
        lcmm = 1;
        MaxTIs = max(MRStruct.Par.nTempIntsPerAngInt);
        if(all(MRStruct.Hdr.Dicom.alICEProgramPara == 0))   % Detect if we are dealing with Lukis original sequence. He always uses as max TI = 3
           MaxTIs = 3;                                        % no matter what real TIs are
        end
        for jk = 2:MaxTIs                
            lcmm = lcm(lcmm,jk);
        end
        if(MRStruct.Par.VTI_Flag)
            MRStruct.Par.vecSize = floor(MRStruct.Par.vecSize  / lcmm) * lcmm; 
        end
        
        
    end
    MRStruct.Par.UsefulADCPtsPerAngInt=MRStruct.Par.vecSize./MRStruct.Par.nTempIntsPerAngInt.*MRStruct.Par.TrajPts;
    MRStruct.Par.ADCdtPerAngInt_ns = MRStruct.Par.Dwelltimes(1).*MRStruct.Par.nTempIntsPerAngInt./MRStruct.Par.TrajPts;       
    % Factor 2 NO LONGER NECESSARY, it's done in read_ascconv_VE11_eh!

    MRStruct.Par.SpatialSpectralEncoding_flag = true;
    MRStruct.Par.SpatialSpectralEncoding_type = 'CRT';

        
        
        
end
