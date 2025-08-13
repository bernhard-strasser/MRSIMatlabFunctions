function [MRStruct,MRStruct_refscan,MRStruct_Noise,AdditionalOut] = io_ReadAverageReshape3DCRTDataOwnRead(file, Settings)
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
    if(~isfield(Settings,'OmitDataSets'))
        Settings.OmitDataSets = {};
    end
    MRStruct_Noise = [];
    % Time total
    TicTotal = tic;
    fprintf('\n\nRead and Reshape data\t\t\t\t...')

    %% Read Data

    MRStructRaw = io_ReadSiemensData(file,Settings);
    
    
%     %% DELETEME DEBUG: I USED WRONG SEQUENCE, SO RENAME FIELD. SHOULD NOT BE DONE NORMALLY 
%     MRStructRaw.Data.RTFEEDBACK = MRStructRaw.Data.PATREFANDIMASCAN; MRStructRaw.Data = rmfield(MRStructRaw.Data,'PATREFANDIMASCAN');
%     MRStructRaw.mdhInfo.RTFEEDBACK = MRStructRaw.mdhInfo.PATREFANDIMASCAN; MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,'PATREFANDIMASCAN');
    

    %% Read Parameters

    fprintf('\tReading Parameters\t\t\t...')
    ticcyReadPars = tic;
    [MRStruct,MRStructRaw] = io_Read3DCRTParsOwnRead(MRStructRaw,'ONLINE'); 
    [MRStruct_refscan,MRStructRaw] = io_Read3DCRTParsOwnRead(MRStructRaw,'PATREFSCAN');
    AdditionalOut.CoilCompScan = io_Read3DCRTParsOwnRead(MRStructRaw,'PATREFANDIMASCAN');
    MRStruct.Par.dicom_flag = false;
    MRStruct_refscan.Par.dicom_flag = false;
    fprintf('\n\t\t\t\t\t\t...\ttook\t%10.6f seconds',toc(ticcyReadPars))


    %% Reshape Data 
    ticcyReadAverageReshape = tic;
    fprintf('\n\tReading, Averaging & Reshaping data\t...\n')
    if(isfield(MRStructRaw.Data,'ONLINE'))
        [MRStruct,MRStructRaw] = Reshape3DCRTDataOwnRead(MRStructRaw,MRStruct,'ONLINE');
    end
    if(nargout > 1 && isfield(MRStructRaw.Data,'PATREFSCAN'))
    	[MRStruct_refscan,MRStructRaw] = Reshape3DCRTDataOwnRead(MRStructRaw,MRStruct_refscan,'PATREFSCAN');
    end
    if(nargout > 2 && isfield(MRStructRaw.Data,'NOISEADJSCAN'))
    	[MRStruct_Noise,MRStructRaw]  = ReshapeNoisePrescan(MRStructRaw);
    end

    if(nargout > 3)
    	[AdditionalOut.CoilCompScan,MRStructRaw] = Reshape3DCRTDataOwnRead(MRStructRaw,AdditionalOut.CoilCompScan,'PATREFANDIMASCAN');
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

	
    if(~isfield(MRStructRaw.mdhInfo,CurDataSet))
        return
    end
    
	InputName = inputname(1);
	if(~isempty(InputName))
	    evalin('caller',['clear ' InputName])
	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape raw Data 1 & Perform averaging                                               %
    % Current Size: {nTotADCs}[nADCPts]                                                    %
    % Target Size:  {nCircles}[nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs x nCha] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Data
    TmpData2 = cell([1 MRStruct.Par.nAngInts]); 
    for CurCrcl = 1:MRStruct.Par.nAngInts
        TmpData2{CurCrcl} = single(zeros([MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC(CurCrcl) MRStruct.Par.nADCsPerAngInt(CurCrcl) MRStruct.Par.total_channel_no_reco MRStruct.Par.nRep])); 
    end
    % Rearrange Data
    for ii = 1:numel(MRStructRaw.Data.(CurDataSet))
            
        CurCrcl = MRStructRaw.mdhInfo.(CurDataSet).Lin(ii);

            
        % Remark: The 3rd circle might have partitions (2,3,...,14) saved in .Seg(ii). However, we don't want to zero-fill the other data here,
        % so our TmpData2 only has 13 partitions. Therefore using .Seg(ii) will result in an error, and we need to use the numbers to 
        % (1,2,...,13), which is saved in CurSegToSave
        CurSegToSave = MRStructRaw.mdhInfo.(CurDataSet).Seg(ii) - (MRStruct.Par.nPartEnc_Meas-size(TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)},2))/2;
        CurData = reshape(MRStructRaw.Data.(CurDataSet){ii},[1 1 1 MRStruct.Par.nPtsPerADC( CurCrcl) 1 MRStruct.Par.total_channel_no_reco]);
        MRStructRaw.Data.(CurDataSet){ii} = [];
        TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)}(MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),CurSegToSave,MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Ida(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Rep(ii)) = ...
        TmpData2{MRStructRaw.mdhInfo.(CurDataSet).Lin(ii)}(MRStructRaw.mdhInfo.(CurDataSet).Idb(ii),CurSegToSave,MRStructRaw.mdhInfo.(CurDataSet).Sli(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Ida(ii),:,MRStructRaw.mdhInfo.(CurDataSet).Rep(ii)) + CurData; 
    end
    MRStructRaw.Data = rmfield(MRStructRaw.Data,CurDataSet);
    for ii = 1:numel(TmpData2)
        TmpData2{ii} = TmpData2{ii}/MRStructRaw.mdhInfo.(CurDataSet).NSet;     % For averaging
    end
    MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,CurDataSet);
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
        TmpData3 = reshape(TmpData2{CurCrcl},[MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC(CurCrcl)*MRStruct.Par.nADCsPerAngInt(CurCrcl) MRStruct.Par.total_channel_no_reco MRStruct.Par.nRep]); 
        TmpData3 = TmpData3(:,:,:,1:MRStruct.Par.UsefulADCPtsPerAngInt(CurCrcl),:);
        TmpData3 = reshape(TmpData3,[MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.TrajPts(CurCrcl) MRStruct.Par.vecSize/MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.total_channel_no_reco MRStruct.Par.nRep] ); 
        TmpData3 = permute(TmpData3,[4 2 3 1 5 6 7]); % Permute from [nCha x nTempIntsPerAngInt x nPart x nSlc x nTrajPoints x vecSize/nTempIntsPerAngInt] to [nTrajPoints x nPart x nSlc x vecSize/nTempIntsPerAngInt x nTempIntsPerAngInt x nCha]
        TmpData3 = reshape(TmpData3,[MRStruct.Par.TrajPts(CurCrcl) 1 MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco MRStruct.Par.nRep]); % Merge vecSize/nTempIntsPerAngInt x nTempIntsPerAngInt to vecSize
        TmpData2{CurCrcl} = [];
        MRStruct.Data{CurCrcl} = TmpData3;
%         for Curkz = 1:MRStruct.Par.nPartEnc
%             MRStruct.Data{Curkz}(:,CurCrcl,:,:,:,:) = TmpData3(:,:,Curkz,:,:,:);
%         end
    end


    MRStruct.Par.DataSize = cellfun(@size,MRStruct.Data,'uni',false);     % Generalize later for more than 1 partition!

    
    MRStruct.Par.dimnames_small_cell = {'nAngInts'};
    
    MRStruct.Par.dimnames_small = {'TrajPts','dummy','Part','Slc','vecSize','cha','Rep'};
    
    
    
    %% Averaging
    MRStruct.RecoPar = MRStruct.Par;
    MRStruct.RecoPar.nAverages = 1;
    MRStruct.RecoPar.AverageDataDuringReshape_flag = true;
    
    
    
end


function [MRStruct,MRStructRaw] = ReshapeNoisePrescan(MRStructRaw)    

	
% InputName = inputname(1);
% if(~isempty(InputName))
%     evalin('caller',['clear ' InputName])
% end

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
    
    if(isfield(MRStruct.Par,'Dwelltimes') && isnan(MRStruct.Par.Dwelltimes))            % Assume standard value if we dont know it. Actually we dont write it to header
        MRStruct.Par.Dwelltimes = 10000;                                                % So we cannot know it currently
    end
    
    MRStruct.Par.DataSize = size(MRStruct.Data);     % Generalize later for more than 1 partition!
    MRStructRaw.mdhInfo = rmfield(MRStructRaw.mdhInfo,CurDataSet);
    MRStructRaw.Data = rmfield(MRStructRaw.Data,CurDataSet);
    
    
    
end




function [MRStruct, MRStructRaw] = io_Read3DCRTParsOwnRead(MRStructRaw,CurDataset)


    if(~isfield(MRStructRaw.mdhInfo,CurDataset) || ~isfield(MRStructRaw.Data,CurDataset))
        MRStruct = struct();
        return;
    end


    % CORRECT MRStructRaw
    
    % Our int32 values that we write into the iceParam gets converted to uint32 by the iceParam. Therefore, all negative values that are written into the header are interpreted wrongly:
    % Negative numbers x<0 are represented as y = 2^16 - x in the bit pattern, and so e.g. -1 is represented as 2^16-1 = 65535, -2 as 65534 etc.
    % In total we can represent 2^16 numbers with 16 bits. We split this representation to representing zero, and then 2^15 positive numbers and (2^15)-1 negative numbers
    % (-2^15-1,...,0,...,2^15). So if any uint representation y of a number x is greater than 2^15, it was in fact negative, and its value was x = y - 2^16
    Logi = MRStructRaw.mdhInfo.(CurDataset).iceParam(1,:) > 2^15;
    MRStructRaw.mdhInfo.(CurDataset).iceParam(1,Logi) = MRStructRaw.mdhInfo.(CurDataset).iceParam(1,Logi) -65536;

    Logi = MRStructRaw.mdhInfo.(CurDataset).iceParam(2,:) > 2^15;
    MRStructRaw.mdhInfo.(CurDataset).iceParam(2,Logi) = MRStructRaw.mdhInfo.(CurDataset).iceParam(2,Logi) -65536;

    Logi = MRStructRaw.mdhInfo.(CurDataset).iceParam(3,:) > 2^15;
    MRStructRaw.mdhInfo.(CurDataset).iceParam(3,Logi) = MRStructRaw.mdhInfo.(CurDataset).iceParam(3,Logi) -65536;
    
    Logi = MRStructRaw.mdhInfo.(CurDataset).iceParam(4,:) > 2^15;
    MRStructRaw.mdhInfo.(CurDataset).iceParam(4,Logi) = MRStructRaw.mdhInfo.(CurDataset).iceParam(4,Logi) -65536;
    
    Logi = MRStructRaw.mdhInfo.(CurDataset).iceParam(6,:) > 2^15;
    MRStructRaw.mdhInfo.(CurDataset).iceParam(6,Logi) = MRStructRaw.mdhInfo.(CurDataset).iceParam(6,Logi) -65536;   

    
    %fn: change hamming addresses to 10,11,12
    HammPartitions = (-1).^MRStructRaw.mdhInfo.(CurDataset).iceParam(10,:) .* (MRStructRaw.mdhInfo.(CurDataset).iceParam(11,:) + MRStructRaw.mdhInfo.(CurDataset).iceParam(12,:)/10000);
           
    if(numel(unique(HammPartitions)) > 1 && 2*MRStructRaw.Par.WipMemBlock_alFree(1) ~= MRStructRaw.Par.nFreqEnc)   % Currently, in Fabians sequence you can do only Hamming
        MRStructRaw.Par.Hamming_flag = true;                                                                       % in all three dimensions. If the number of nFreqEnc 
    else                                                                                                           % is not what the sequence wrote as being the number of
        MRStructRaw.Par.Hamming_flag = false;                                                                      % circles, Hamming must have been enabled
    end
    
    
    if(~strcmpi(CurDataset,'NoiseAdjScan') && MRStructRaw.Par.nPartEnc > 1 && (strcmpi(MRStructRaw.Par.BaselineVersion,'vb') || all(MRStructRaw.mdhInfo.(CurDataset).Seg == 1) ))
        MRStructRaw = CorrMdhForOldADC(MRStructRaw,CurDataset);
    end
    
    
    
    
    % Copy all fields
    for CurField = transpose(fieldnames(MRStructRaw))
        if(strcmpi(CurField{1},'Data') || strcmpi(CurField{1},'mdhInfo'))
            continue;
        end
        MRStruct.(CurField{1}) = MRStructRaw.(CurField{1});
    end
    
    MRStruct.Par.nAngInts = MRStructRaw.mdhInfo.(CurDataset).NLin;
    
    % Find first occurence of circle CurCrcl. Take the info for that circle. All the other circles have the same entries and are not needed
    for CurCrcl = 1:MRStruct.Par.nAngInts   
        FirstCircleMeas = find(MRStructRaw.mdhInfo.(CurDataset).Lin == CurCrcl,1); 
        MRStruct.Par.FirstCirclekSpacePoint(1,CurCrcl) = (MRStructRaw.mdhInfo.(CurDataset).iceParam(1,FirstCircleMeas) + MRStructRaw.mdhInfo.(CurDataset).iceParam(2,FirstCircleMeas)/10000);
        MRStruct.Par.FirstCirclekSpacePoint(2,CurCrcl) = (MRStructRaw.mdhInfo.(CurDataset).iceParam(3,FirstCircleMeas) + MRStructRaw.mdhInfo.(CurDataset).iceParam(4,FirstCircleMeas)/10000);
    end
    

    MRStruct.Par.nRep = max(MRStructRaw.mdhInfo.(CurDataset).Rep); 

    
    MRStruct.Par.ReadoutOSFactor = MRStructRaw.Hdr.Dicom.flReadoutOSFactor;
%     MRStruct.Par.Dwelltimes = MRStruct.Par.Dwelltimes/MRStruct.Par.ReadoutOSFactor;       % Alrdy corrected in read_ascconv_VE11_eh
    
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

    % For getting vecSize, it depends on whether we read in ONLINE or refscan data
    if(strcmpi(CurDataset,'PATREFANDIMASCAN'))
        if(any(MRStructRaw.Hdr.Dicom.alICEProgramPara(18:21)))
            MRStruct.Par.nFreqEnc = MRStructRaw.Hdr.Dicom.alICEProgramPara(20);
            MRStruct.Par.nPhasEnc = MRStruct.Par.nFreqEnc;
            MRStruct.Par.nPartEnc = MRStructRaw.Hdr.Dicom.alICEProgramPara(21); 
            MRStruct.Par.Dwelltimes = MRStructRaw.Hdr.Dicom.alICEProgramPara(18)*1E3; 
            MRStruct.Par.vecSize = MRStruct.Hdr.Dicom.alICEProgramPara(19);  
        else
            MRStruct.Par.vecSize = 5;  
            if(isfield_recursive(MRStructRaw,'mdhInfo.ONLINE.Idb'))
                MRStruct.Par.Dwelltimes = MRStruct.Par.Dwelltimes/max(MRStructRaw.mdhInfo.ONLINE.Idb);
            end
        end
    elseif(strcmpi(CurDataset,'ONLINE'))
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
    
    % For Lukas's sequence, he writes the TrajPts into Idc. Last if condition: In Lukis newer sequence, the IceParam(6,:) stores the partition number 
    if(any(MRStruct.Par.TrajPts == 0) || all(MRStruct.Par.TrajPts == MRStruct.Par.ReadoutOSFactor) || all(MRStruct.Par.TrajPts + 1 == MRStruct.Par.nPartEncsPerAngInt))
        MRStruct.Par.TrajPts = repmat(MRStruct.Par.ReadoutOSFactor*(MRStructRaw.mdhInfo.(CurDataset).Idc(1)-1),[1 MRStruct.Par.nAngInts]); 
    end

    % Define which angular interleaves are measured for each kz (matrix of nPartEncmax x nCircles)
%     MRStruct.Par.nPartEncsPerAngInt = [9 9 9 7 7 7 7 5 5 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1]; % For simulating 3D
    dummy1 = max(MRStruct.Par.nPartEncsPerAngInt);
    dummy2 = ceil(dummy1/2);
    MRStruct.Par.kzPositions = zeros([1 dummy1]);
    for kz = 1:max(MRStruct.Par.nPartEncsPerAngInt)
        MRStruct.Par.AngIntsPerPartEnc(kz,:) = MRStruct.Par.nPartEncsPerAngInt >= dummy1-(     (kz - (kz>dummy2)*2*mod(kz,dummy2))        -1)*2;
        
        % Restrict somehow to only define this in certain cases?
        %fn: change hamming addresses to 10,11,12
        MRStruct.Par.kzPositions(kz) = (-1).^max(MRStructRaw.mdhInfo.(CurDataset).iceParam(10,MRStructRaw.mdhInfo.(CurDataset).Seg == kz)) .* (max( MRStructRaw.mdhInfo.(CurDataset).iceParam(11,MRStructRaw.mdhInfo.(CurDataset).Seg == kz)) + max(MRStructRaw.mdhInfo.(CurDataset).iceParam(12,MRStructRaw.mdhInfo.(CurDataset).Seg == kz)/10000));

        
    end

    % Hack -- Delete me later
    if(MRStruct.Par.TrajPts == 0)
        MRStruct.Par.TrajPts = [4 10 18 24 30 36 45 48 60 60 72 80 80 90 120 120 120 120 120 144 144 144 144 180 180 180 180 180 180 216 216 216] * 2;
    end

    MRStruct.Par.nPartEnc_Meas = max(MRStruct.Par.nPartEncsPerAngInt);  % The measured ones can be different than the one written in the header

    MRStruct.Par.RewPts = 0;
    MRStruct.Par.TrajTotPts = MRStruct.Par.TrajPts + MRStruct.Par.RewPts;

    
    % Three-D Hamming filter
    % In very rare cases this could be true although there was no 3D Hamming measured...
    MRStruct.Par.ThreeDHamming_flag = MRStructRaw.Par.Hamming_flag && (MRStruct.Par.nPartEnc_Meas > MRStruct.Par.nPartEnc);

    MRStruct.Par.VTI_Flag = ~all(MRStruct.Par.nTempIntsPerAngInt == MRStruct.Par.nTempIntsPerAngInt(1));

    %coil compression
    MRStruct.Par.coilcompression_to_numb_coils=MRStructRaw.Hdr.Dicom.alICEProgramPara(10);
    % MRStruct.Par.coilcompression_to_numb_coils = 8;
    if (MRStruct.Par.coilcompression_to_numb_coils<1 || MRStruct.Par.coilcompression_to_numb_coils>MRStruct.Par.total_channel_no_reco)
        MRStruct.Par.coilcompression_to_numb_coils=MRStruct.Par.total_channel_no_reco;
    end


    MRStruct.Par.UsefulADCPtsPerAngInt=MRStruct.Par.vecSize./MRStruct.Par.nTempIntsPerAngInt.*MRStruct.Par.TrajPts;
    
    if(isfield(MRStruct.Par,'WipMemBlock_alFree') && numel(MRStruct.Par.WipMemBlock_alFree) > 58)
        fprintf('\n\n\n\nWARNING: New Dwelltime-calculation based on wipmemblock!\nOld Dwelltime: %f.',MRStruct.Par.Dwelltimes)
%         MRStruct.Par.Dwelltimes = MRStruct.Par.WipMemBlock_alFree(59)/MRStruct.Par.ReadoutOSFactor*MRStruct.Par.TrajPts(1)/MRStruct.Par.nTempIntsPerAngInt(1);   % ADC-dt times TrajPts divided by TempInt
        MRStruct.Par.Dwelltimes = round(MRStruct.Par.Dwelltimes/10^4)*10^4;
        fprintf('\nNew  Dwelltime: %f.\nIf there is a discrepancy, the reco might be wrong.\n\n\n\n',MRStruct.Par.Dwelltimes)
    end
    
    MRStruct.Par.ADCdtPerAngInt_ns = MRStruct.Par.Dwelltimes(1).*MRStruct.Par.nTempIntsPerAngInt./MRStruct.Par.TrajPts;       
    % Factor 2 NO LONGER NECESSARY, it's done in read_ascconv_VE11_eh!

    MRStruct.Par.SpatialSpectralEncoding_flag = true;
    MRStruct.Par.SpatialSpectralEncoding_type = 'CRT';

        
        
        
end



%%
function MRStruct = CorrMdhForOldADC(MRStruct,CurSet)

		
	InputName = inputname(1);
	if(~isempty(InputName))
	    evalin('caller',['clear ' InputName])
	end

    Settings = struct;
    if(strcmpi(MRStruct.Par.AssumedSequence,'ViennaCRT'))
        Tmpp = numel(MRStruct.mdhInfo.(CurSet).Col)/MRStruct.mdhInfo.(CurSet).NIda;
        if(abs(Tmpp - round(Tmpp)) > 0)
            MissingPts = round(abs(Tmpp - round(Tmpp)) * MRStruct.mdhInfo.(CurSet).NIda);
            for CurField = {'Col','Lin','Ave','Sli','Par','Eco','Phs','Rep','Set','Seg','Idb','Idc','Idd','Ide','iceParam'}
                MRStruct.mdhInfo.(CurSet).(CurField{1}) = cat(2,MRStruct.mdhInfo.(CurSet).(CurField{1}),MRStruct.mdhInfo.(CurSet).(CurField{1})(:,end-MissingPts+1:end));
            end
            AppendIda = MRStruct.mdhInfo.(CurSet).Ida(end)+1 : (MRStruct.mdhInfo.(CurSet).Ida(end)-MRStruct.mdhInfo.(CurSet).Ida(end-1)) : (MRStruct.mdhInfo.(CurSet).Ida(end)+MissingPts);
            MRStruct.mdhInfo.(CurSet).Ida = cat(2,MRStruct.mdhInfo.(CurSet).Ida,AppendIda);
        end
    end


    Crcls = (MRStruct.mdhInfo.(CurSet).iceParam(1,:) + MRStruct.mdhInfo.(CurSet).iceParam(2,:)/10000);
    NCircl = numel(unique( Crcls ));
    
    if(MRStruct.Par.Hamming_flag)
        %fn: change hamming addresses to 10,11,12
        HammPartitions = (-1).^MRStruct.mdhInfo.(CurSet).iceParam(10,:) .* (MRStruct.mdhInfo.(CurSet).iceParam(11,:) + MRStruct.mdhInfo.(CurSet).iceParam(12,:)/10000);
        NPart = numel(unique(HammPartitions)); 
    else                                                                                                    
        NPart = MRStruct.Par.nPartEnc;
    end
    NPartMax = floor(NPart/2);

    
    
    % Calc CirclRadius and Partitions and kDist
    CirclRadius_dummy = transpose((0:(NCircl-1)) +0.5);

    Part = [];
    CirclRadius = [];
    for ii = -NPartMax:NPartMax
        CirclesPerPart(ii+NPartMax+1) = max(2,ceil(sqrt(NCircl^2-ii^2*NCircl^2/NPartMax^2)));

        Part = cat(1,Part,repmat(ii,[CirclesPerPart(ii+NPartMax+1) 1]));
        CirclRadius = cat(1,CirclRadius,CirclRadius_dummy(1:CirclesPerPart(ii+NPartMax+1)));
    end

%     kxy = CirclRadius * Deltakxy; kz = Part * Deltakz;
%     kDist = sqrt(kxy.^2 + kz.^2); 
% 
% 
%     % Sort like Luki does
% 
%     pattern = 1:numel(Part);
% 
%     linind = 0;
%     for ii = 1:numel(Part)
%         for jj = ii:numel(Part)
%             
%             if(ii == 14 && jj == 170)
%                 bla = 1;
%             end
% 
%             if(kDist(jj) < kDist(ii))
%                 linind = linind + 1;
%                 ii_matlab(linind) = ii; jj_matlab(linind) = jj;
%                 jTemp = pattern(ii); 
%                 pattern(ii) = pattern(jj);
%                 pattern(jj) = jTemp;
%                 kTemp = kDist(ii);
%                 kDist(ii) = kDist(jj);
%                 kDist(jj) = kTemp;
%             end
%         end
%     end
% 
% 
%     % Create Vector containing indices, because every ADC, TempInt, Average (Set) and repetition (Rep) has the same value
%     Indie = cumsum(MRStruct.mdhInfo.(CurSet).Ida == 1 & MRStruct.mdhInfo.(CurSet).Idb == 1 & MRStruct.mdhInfo.(CurSet).Set == 1 & MRStruct.mdhInfo.(CurSet).Rep == 1);
% 
%     % Replicate the pattern, because each ADC etc. has the same pattern value (as long as nothing else changes)
%     pattern = pattern(Indie);


    % The calculation with computing the radius does sometimes not work, because of digital errors (in the order of 1E-16), where the compiler in IDEA is saying kDist(jj) < kDist(ii), while matlab says
    % kDist(jj) > kDist(ii), when they are actually identical up to numerical errors. So instead I calculate the order of circles and partitions differently, by knowing the loops per partition, and knowing that
    % Luki always counts from the north pole down to the south pole with the linear indices (so northpole partition + loop 1 is linear index 1, then northpole partition + loop 2 is 2, then northpole-1 partition
    % loop 1 is 3 etc.).
    if(MRStruct.Par.Hamming_flag)
        MRStruct.mdhInfo.(CurSet).Seg = cat(2,1,cumsum(diff(HammPartitions) > 0) + 1);
        MRStruct.mdhInfo.(CurSet).Lin = ReplaceValuesInNumArray(Crcls,unique(Crcls),1:numel(unique(Crcls)));
    else
        CurInd = 1; 
        for ii = 1:numel(CirclesPerPart)
            LinNew(CurInd:CurInd+CirclesPerPart(ii)-1) = 1:CirclesPerPart(ii); 
            SegNew(CurInd:CurInd+CirclesPerPart(ii)-1) = ii; 
            CurInd = CurInd + CirclesPerPart(ii); 
        end
        LinNew = LinNew(MRStruct.mdhInfo.(CurSet).Lin);
        MRStruct.mdhInfo.(CurSet).Seg = SegNew(MRStruct.mdhInfo.(CurSet).Lin);
        MRStruct.mdhInfo.(CurSet).Lin = LinNew;
    end
    
    % Overwrite Seg and Lin. Seg should contain the current partition, Lin, the current circle
%     MRStruct.mdhInfo.(CurSet).Lin = ceil(CirclRadius(pattern));
%     MRStruct.mdhInfo.(CurSet).Seg = floor(Part(pattern)) + NPartMax+1;
    MRStruct.mdhInfo.(CurSet).NLin = max(MRStruct.mdhInfo.(CurSet).Lin);
    MRStruct.mdhInfo.(CurSet).NSeg = max(MRStruct.mdhInfo.(CurSet).Seg);

    MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);

end

