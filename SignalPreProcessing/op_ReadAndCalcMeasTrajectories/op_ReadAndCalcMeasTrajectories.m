function [kSpaceTrajectory,kSpaceTrajectoryAllFIDPts] = op_ReadAndCalcMeasTrajectories(InputFile,Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [Output, AdditionalOut] = op_ReconstructNonCartMRData(Output,AdditionalIn,Settings)
%
% Input: 
% -         ?                     ...     
% -         ?                     ...     
% -         ?             ...     
%
% Output:
% -         ?                      ...     
% -         ?                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations


if(~exist('AdditionalIn','var'))
    AdditionalIn = struct();
end
if(~exist('Settings','var'))
    Settings = struct();
end
if(~isfield(Settings,'Plot'))
    Settings.Plot = struct;
end
if(~isfield(Settings.Plot,'PlotAveragedTraj'))
    Settings.Plot.PlotAveragedTraj = false;
end

if(~isfield(Settings,'m_lADCTimeShift_us'))
    Settings.m_lADCTimeShift_us = 1; 
end
if(~isfield(Settings,'RemoveFirstFractionOfTrajPtsFromBeginning'))  % Because the first circumnaviagion is often wrong. Should be between 0 and 1, e.g. 0.2.
    Settings.RemoveFirstFractionOfTrajPtsFromBeginning = 0; 
end                       
                      
% Output = supp_FixPars(Output);  % To hard-code/hack parameters for special cases, or to make Parameters consistent between different read-in-methods.


tic
fprintf('\n\nReading and calculating trajectory took \t\t...')




%%

Tmpp = read_ascconv(InputFile);

if(strcmpi(Tmpp.AssumedSequence,'ViennaCRT'))
    MeasData = io_ReadAverageReshape3DCRTDataOwnRead_TrajMeas(InputFile);
else
    MeasData = io_ReadAverageReshapeGenericData(InputFile);
end


%%



if(strcmpi(Tmpp.AssumedSequence,'ViennaCRT'))
    Tmp = MeasData; Tmp.Data = [];

    for ii = 1:numel(MeasData.Data)

        TrajPts = Tmp.Par.TrajPts(ii);
        StartingPointAfterLaunchTrack = round(TrajPts*Settings.RemoveFirstFractionOfTrajPtsFromBeginning);
        LaunchPts = Tmp.Par.LaunchTrackADCPts(ii) + StartingPointAfterLaunchTrack;
        nTempInt = Tmp.Par.nTempIntsPerAngInt(ii);
        vecSize = Tmp.Par.vecSize/nTempInt;

        Tmp.Data = permute(MeasData.Data{ii},[1 5 3 8 2 4 6 7]);

        Sizz = size(Tmp.Data); Tmp.Data = reshape(Tmp.Data,[Sizz(1)*Sizz(2) 1 Sizz(3:end)]); Tmp2 = op_CalcTraj(Tmp);

%             Tmp2 = op_CalcTraj(Tmp);
%             Sizz = size(Tmp2.x); if(numel(Sizz) < 3); Sizz(3) = 1; end
%             Tmp2.x = reshape(Tmp2.x,[Sizz(1)*Sizz(2) Sizz(3:end)]); Tmp2.y = reshape(Tmp2.y,[Sizz(1)*Sizz(2) Sizz(3:end)]);


        Tmp2.x = Tmp2.x(1:end-(2*TrajPts-LaunchPts),:,:,:);
        Tmp2.y = Tmp2.y(1:end-(2*TrajPts-LaunchPts),:,:,:);


        Tmp2.x= reshape(Tmp2.x(LaunchPts+1:end,:,:,:),[TrajPts vecSize-2 1 nTempInt]);
        Tmp2.y= reshape(Tmp2.y(LaunchPts+1:end,:,:,:),[TrajPts vecSize-2 1 nTempInt]);

        k.x{ii} = Tmp2.x; k.y{ii} = Tmp2.y;

        k.StartingPointAfterLaunchTrack{ii} = StartingPointAfterLaunchTrack+1;      % +1 bc that's the actual start where we start with the trajectory


    end

else

    MeasData.Data = permute(MeasData.Data,[1 12 13 14 3 2 4 5 6 7 8 9 10 11]);
    DataSize = size(MeasData.Data); DataSizeOrig = DataSize;
    MeasData.Data = reshape(MeasData.Data,[DataSize(1)*DataSize(2) DataSize(3) DataSize(4) DataSize(5)]);
    MeasData.Data = MeasData.Data(:,:,2:end,:);
    DataSize = size(MeasData.Data);

    TrajPts = MeasData.Par.IceParam(15,1) * MeasData.Par.ReadoutOSFactor;
    nTempInt = MeasData.Par.IceParam(20,1);
    vecSize = DataSize(1)/TrajPts;
%         Crcls = numel(MeasData.Data)/TrajPts/nTempInt/vecSize/8;
    Crcls = DataSize(4);
    ADC_dt = MeasData.Par.Dwelltimes(1)*nTempInt*1E-9/TrajPts;
    LaunchPts_FromHeader = round(MeasData.Par.IceParam(18,1:prod(DataSizeOrig(2:3))*8:end)*10/(ADC_dt*1E6));
%         LaunchPts = LaunchPts + floor(10/(ADC_dt*1E6)*m_lADCTimeShift_us) - 1 + 42;
    StartingPointAfterLaunchTrack = floor(10/(ADC_dt*1E6)*Settings.m_lADCTimeShift_us) - 1 + round(Settings.RemoveFirstFractionOfTrajPtsFromBeginning*TrajPts);
    LaunchPts = LaunchPts_FromHeader + StartingPointAfterLaunchTrack;
    MeasData.Data = reshape(MeasData.Data,[TrajPts vecSize nTempInt 8 Crcls]);
    DataSize = size(MeasData.Data);


    for ii = 1:DataSize(5)
        Tmp = MeasData; Tmp.Data = Tmp.Data(:,:,:,:,ii,:,:);

        Sizz = size(Tmp.Data); Tmp.Data = reshape(Tmp.Data,[Sizz(1)*Sizz(2) 1 Sizz(3:end)]); Tmp2 = op_CalcTraj(Tmp);


%             Tmp2 = op_CalcTraj(Tmp);
%             Sizz = size(Tmp2.x); if(numel(Sizz) < 3); Sizz(3) = 1; end
%             Tmp2.x = reshape(Tmp2.x,[Sizz(1)*Sizz(2) Sizz(3:end)]); Tmp2.y = reshape(Tmp2.y,[Sizz(1)*Sizz(2) Sizz(3:end)]);

        Tmp2.x = Tmp2.x(1:end-(TrajPts-LaunchPts(ii)),:,:,:);
        Tmp2.y = Tmp2.y(1:end-(TrajPts-LaunchPts(ii)),:,:,:);

        % Reshape Data
        k.x{ii} = reshape(Tmp2.x(LaunchPts(ii)+1:end,:,:,:),[TrajPts vecSize-1 1 nTempInt]);
        k.y{ii} = reshape(Tmp2.y(LaunchPts(ii)+1:end,:,:,:),[TrajPts vecSize-1 1 nTempInt]);

        k.StartingPointAfterLaunchTrack{ii} = StartingPointAfterLaunchTrack+1;

    end



end

%     % k-Space Distance
%     for ii = 1:numel(k.x) 
%         k.DistToCtr{ii} = sqrt(k.x{ii}.^2 + k.y{ii}.^2); 
%     end

kAvg2To20 = k;
for ii = 1:numel(k.x) 
    kAvg2To20.x{ii} = mean(k.x{ii}(:,2:20,:,:),2);
    kAvg2To20.y{ii} = mean(k.y{ii}(:,2:20,:,:),2);

end

Trajj = cell([1 numel(kAvg2To20.y)]);
TrajjAllFIDPts = Trajj;
for ii = 1:numel(kAvg2To20.y) 
    Sizz = size(kAvg2To20.y{ii}); Sizz = cat(2,Sizz,ones([1 4-numel(Sizz)]));
    Trajj{ii} = cat(1,-reshape(kAvg2To20.y{ii}(:,:,1,1),[1 Sizz(1) Sizz(2) 1]),-reshape(kAvg2To20.x{ii}(:,:,1,1),[1 Sizz(1) Sizz(2) 1])); 
    Trajj{ii} = Trajj{ii}/(2*MeasData.Par.nFreqEnc/2*2*pi/MeasData.Par.FoV_Read); % Originall it is in /mm, need to normalize it. 2: To get from -0.5 to 0.5. 
 
    Sizz2 = size(k.y{ii}); Sizz2 = cat(2,Sizz2,ones([1 4-numel(Sizz2)]));
    TrajjAllFIDPts{ii} = cat(1,-reshape(k.y{ii},[1 Sizz2]),-reshape(k.x{ii},[1 Sizz2])); 
    TrajjAllFIDPts{ii} = TrajjAllFIDPts{ii}/(2*MeasData.Par.nFreqEnc/2*2*pi/MeasData.Par.FoV_Read); % Originall it is in /mm, need to normalize it. 2: To get from -0.5 to 0.5.
    
    
end
clear kSpaceTrajectory; 
kSpaceTrajectory.StartingPointAfterLaunchTrack = kAvg2To20.StartingPointAfterLaunchTrack;
kSpaceTrajectory.GM = Trajj;

    
kSpaceTrajectoryAllFIDPts.StartingPointAfterLaunchTrack = kAvg2To20.StartingPointAfterLaunchTrack;
kSpaceTrajectoryAllFIDPts.GM = TrajjAllFIDPts;
       
    




%% Postparations

Output = supp_UpdateRecoSteps(Output,Settings);

fprintf('\n\t\t\t\t...\ttook\t%10.6f seconds',toc)


