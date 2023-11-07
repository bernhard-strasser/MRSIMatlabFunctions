function [MRStruct] = op_Reshape3DCRTData(MRStruct, mdhInfo, Settings)
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

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);

% Time Reshaping
tic
fprintf('\n\nReshaping data\t\t\t...')


%% Reshape raw Data 1 & Perform averaging

% Initialize Data
TmpData2 = cell([1 MRStruct.Par.nAngInts]); 
for CurCrcl = 1:MRStruct.Par.nAngInts
    TmpData2{CurCrcl} = zeros([MRStruct.Par.total_channel_no_reco MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC MRStruct.Par.nADCsPerAngInt(CurCrcl) ]); 
end

% Get data. Perform memory-efficient averaging
% TODO: Create NoiseData From Averages
% Current Size {nCircles}[nCha x nTempIntsPerAngInt x nPart x nSlc x nADCPts x nADCs]
for ii = 1:size(MRStruct.Data,3)
    CurData = reshape(transpose(MRStruct.Data(:,:,ii)),size( TmpData2{MRStruct.MapVBVDHeader.Ide(ii)}(:,1,1,1,:,1)));
    TmpData2{MRStruct.MapVBVDHeader.Ide(ii)}(:,MRStruct.MapVBVDHeader.Idb(ii),MRStruct.MapVBVDHeader.Seg(ii),MRStruct.MapVBVDHeader.Sli(ii),:,MRStruct.MapVBVDHeader.Ida(ii)) = TmpData2{MRStruct.MapVBVDHeader.Ide(ii)}(:,MRStruct.MapVBVDHeader.Idb(ii),MRStruct.MapVBVDHeader.Seg(ii),MRStruct.MapVBVDHeader.Sli(ii),:,MRStruct.MapVBVDHeader.Ida(ii)) + CurData; 
end
MRStruct.Data = [];
MRStruct.MapVBVDHeader = [];


%% Reshape Data 2
% Reshape Data to Output Size {MRStruct.Par.nAngInts}[nTrajPoints x nPart x nSlc x vecSize x nCha]

for Curkz = 1:MRStruct.Par.nPartEnc
    MRStruct.Data{Curkz} = zeros([MRStruct.Par.TrajPts MRStruct.Par.nAngInts MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]);
end
% for CurCrcl = 1:MRStruct.Par.nAngInts
%     MRStruct.Data{CurCrcl} = zeros([MRStruct.Par.TrajPts MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]);
% end
for CurCrcl = 1:MRStruct.Par.nAngInts
    TmpData3 = reshape(TmpData2{CurCrcl},[MRStruct.Par.total_channel_no_reco MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.nPtsPerADC*MRStruct.Par.nADCsPerAngInt(CurCrcl) ]); 
    TmpData3 = TmpData3(:,:,:,:,1:MRStruct.Par.UsefulADCPtsPerAngInt(CurCrcl));
    TmpData3 = reshape(TmpData3,[MRStruct.Par.total_channel_no_reco MRStruct.Par.nTempIntsPerAngInt(CurCrcl) MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.TrajPts MRStruct.Par.vecSize/MRStruct.Par.nTempIntsPerAngInt(CurCrcl)]); 
    TmpData3 = permute(TmpData3,[5 3 4 2 6 1]); % New size : [nTrajPoints x nPart x nSlc x vecSize x nTempIntsPerAngInt x nCha]
    TmpData3 = reshape(TmpData3,[MRStruct.Par.TrajPts MRStruct.Par.nPartEncsPerAngInt(CurCrcl) MRStruct.Par.nSLC MRStruct.Par.vecSize MRStruct.Par.total_channel_no_reco]);
%     MRStruct.Data{CurCrcl} = TmpData3;
    for Curkz = 1:MRStruct.Par.nPartEnc
        MRStruct.Data{Curkz}(:,CurCrcl,:,:,:,:) = TmpData3(:,:,Curkz,:,:,:);
    end
end


% Reshape to [nTrajPoints x nAngInt x nPart x nSlc x nTempInt*vecSize x nCha]
Size = size(MRStruct.Data); Size = cat(2,Size,ones([1 7-numel(Size)]));
MRStruct.Data = reshape(MRStruct.Data,[Size(1:4) prod(Size(5:6)) Size(7:end)]);
if(isfield(MRStruct,'NoiseData') && ~isempty(MRStruct.NoiseData))
    MRStruct.NoiseData = reshape(MRStruct.NoiseData,[Size(1:4) prod(Size(5:6)) Size(7:end)]);
end

MRStruct.Par.DataSize = size(MRStruct.Data);


%% Postparations

% Output.Data = MRStruct.Data;
% if(isfield(MRStruct,'NoiseData') && ~isempty(MRStruct.NoiseData))
%     Output.NoiseData = MRStruct.NoiseData;
% end
% clear MRStruct
% Output.MRStruct.Par = MRStruct.Par;

%% The End

fprintf('\n\t\t\t\t...\ttook\t%10.6f seconds',toc)


