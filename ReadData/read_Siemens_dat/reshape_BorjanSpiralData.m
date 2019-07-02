function [kSpace,Par] = reshape_BorjanSpiralData(kSpace, mdhInfo, Par, Settings)

%% Prep
tic
fprintf('\n\nReshaping data\t\t\t...')

if(~exist('Settings','var') || ~isfield(Settings,'ProduceNoiseThroughAvgs_flag'))
    Settings.ProduceNoiseThroughAvgs_flag = true;
end

%% 
kSpace = reshape_Siemens_dat(kSpace,mdhInfo,[11 9 2 1 7 4 5],true); %,[1 7 4 5 11 9 2]); 
kSpace.ONLINE{1} = single(kSpace.ONLINE{1});
clear mdhInfo

% Reshape to [nTempInt x nAngInt x samples x nADCs x nCha x nAvg x nPart x nSlc]
CurSize = size(kSpace.ONLINE{1}); CurSize = cat(2,CurSize,ones([1 8-numel(CurSize)]));
kSpace.ONLINE{1} = reshape(kSpace.ONLINE{1},[Par.nTempInt CurSize(1)/Par.nTempInt CurSize(2:end)]);
% kSpace.ONLINE{1} = reshape(kSpace.ONLINE{1},[CurSize(1)/Par.nTempInt Par.nTempInt CurSize(2:5) prod(CurSize(6:7))]);


% Reshape to [nTempInt x nAngInt x samples*nADCs x nCha x nAvg x nPart x nSlc]
% Cut those ADC-points that are noise only away
ADC_Points = Par.vecSize * Par.TrajTotPts/Par.nTempInt;        % 
Size_D2 = size(kSpace.ONLINE{1});
kSpace.ONLINE{1} = reshape(kSpace.ONLINE{1},[Size_D2(1:2) prod(Size_D2(3:4)) Size_D2(5:end)]);


% Cut together the useful data from the different temp interleaves
% (for the first TI the points 1:ADC_Points are useful, for the second TI the points Par.TrajTotPts/Par.nTempInt+1:ADC_Points+Par.TrajTotPts/Par.nTempInt are useful,
% for the third TI the points 2*Par.TrajTotPts/Par.nTempInt+1:ADC_Points+2*Par.TrajTotPts/Par.nTempInt are useful etc.
PtNumber = size(kSpace.ONLINE{1},3) - ADC_Points;
LastPts = kSpace.ONLINE{1}(:,:,end-PtNumber+1:end,:,:,:,:);
kSpace.ONLINE{1} = kSpace.ONLINE{1}(:,:,1:ADC_Points,:,:,:,:);
for curTempInt = 1:Par.nTempInt
    StartPt = (curTempInt-1)/Par.nTempInt*Par.TrajTotPts;        
%     EndPt = ADC_Points - StartPt;
    kSpace.ONLINE{1}(curTempInt,:,:,:,:,:,:) = cat(3,kSpace.ONLINE{1}(curTempInt,:,StartPt+1:ADC_Points,:,:,:,:), LastPts(curTempInt,:,1:StartPt,:,:,:,:));
end
clear LastPts

% Reshape to [nTempInt x nAngInt x nTrajPoints x vecSize x nCha x nAvg x nPart x nSlc]
% Cut the rewinder away, and save the final size
kSpace.ONLINE{1} = reshape(kSpace.ONLINE{1},[Size_D2(1:2) Par.TrajTotPts Par.vecSize/Par.nTempInt Size_D2(5:end)]);
kSpace.ONLINE{1} = kSpace.ONLINE{1}(:,:,1:Par.TrajPts,:,:,:,:,:);


% % DEBUG: ONLY USE FIRST FID POINT
% kSpace.ONLINE{1} = kSpace.ONLINE{1}(:,1,:,1,1,1);
% Par.vecSize = 1; Par.nTempInt = 1; Par.total_channel_no_measured = 1; Par.vecSize = 1;


% Define Noise
if(~isfield(kSpace,'NOISEADJSCAN') && Settings.ProduceNoiseThroughAvgs_flag && size(kSpace.ONLINE{1},6) > 1)
    kSpace.NOISEADJSCAN{1} = kSpace.ONLINE{1}(:,:,:,:,:,2,:,:) - kSpace.ONLINE{1}(:,:,:,:,:,1,:,:);
    kSpace.NOISEADJSCAN{1} = kSpace.NOISEADJSCAN{1} * sqrt(size(kSpace.ONLINE{1},6)/2); % Scale the noise to be the same as the noise in the actual data is
else
    kSpace.NOISEADJSCAN{1} = [];
end
    
% Add averages together
kSpace.ONLINE{1} = squeeze_single_dim(sum(kSpace.ONLINE{1},6),6);  % squeeze causes non-problematic error in case there is no part/slc
% kSpace.ONLINE{1} = kSpace.ONLINE{1}(:,:,:,:,:,1,:,:);

% Permute from [nTempInt x nAngInt x nTrajPoints x vecSize x nCha x nPart x nSlc] 
%           to [nTrajPoints x nAngInt x nPart x nSlc x nTempInt x vecSize x nCha]
kSpace.ONLINE{1} = permute(kSpace.ONLINE{1},[3 2 6 7 1 4 5]);
if(isfield(kSpace,'NOISEADJSCAN') && ~isempty(kSpace.NOISEADJSCAN{1}))
    kSpace.NOISEADJSCAN{1} = permute(kSpace.NOISEADJSCAN{1},[3 2 6 7 1 4 5]);
end

% Reshape to [nTrajPoints x nAngInt x nPart x nSlc x nTempInt*vecSize x nCha]
Size = size(kSpace.ONLINE{1}); Size = cat(2,Size,ones([1 7-numel(Size)]));
kSpace.ONLINE{1} = reshape(kSpace.ONLINE{1},[Size(1:4) prod(Size(5:6)) Size(7:end)]);
if(isfield(kSpace,'NOISEADJSCAN') && ~isempty(kSpace.NOISEADJSCAN{1}))
    kSpace.NOISEADJSCAN{1} = reshape(kSpace.NOISEADJSCAN{1},[Size(1:4) prod(Size(5:6)) Size(7:end)]);
end



%% The End

fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)


