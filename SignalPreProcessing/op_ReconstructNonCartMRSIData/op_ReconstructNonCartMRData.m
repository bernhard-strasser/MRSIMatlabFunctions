function [Output, AdditionalOut] = op_ReconstructNonCartMRData(Output,AdditionalIn,Settings)
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
if(~isfield(Settings,'Phaseroll_flag'))
    Settings.Phaseroll_flag = true;    
end
if(~isfield(Settings,'DensComp_flag'))
    Settings.DensComp_flag = true;    
end
if(~isfield_recursive(Settings,'DensComp.Method'))
   Settings.DensComp.Method = 'AntoinesVoronoi'; 
%    Settings.DensComp.Method = 'SpiralHoge1997AbruptChanges'; 
end
if(~isfield_recursive(Settings,'DensComp.AutoScale_flag'))
    Settings.DensComp.AutoScale_flag = false;    
end
if(~isfield(Settings.DensComp,'Normalize_flag'))
    Settings.DensComp.Normalize_flag = false;    
end
if(~isfield(Settings,'ConjInkSpace_flag'))
   Settings.ConjInkSpace_flag = false;    
end
if(~isfield(Settings,'ConjIniSpace_flag'))
   Settings.ConjIniSpace_flag = true;    
end
if(~isfield(Settings,'Correct4SpatialB0_flag'))
   Settings.Correct4SpatialB0_flag = false;    
end
if(~isfield(Settings,'CircularSFTFoV_flag'))
   Settings.CircularSFTFoV_flag = false;    
end
if(~isfield(Settings,'PerformZFFT_flag'))
    Settings.PerformZFFT_flag = true;
end
if(~isfield(Settings,'ExtraFoVShiftInPatCoordSys'))
    Settings.ExtraFoVShiftInPatCoordSys = [0 0 0];
end
if(~isfield(Settings,'FlipDim_flag'))
    Settings.FlipDim_flag = false;
end
if(~isfield(Settings,'FlipDim'))
    Settings.FlipDim = 1;
end
if(exist('AdditionalIn','var') && isfield(AdditionalIn,'B0'))
    B0 = AdditionalIn.B0;
end
if(~isfield(Output,'RecoPar'))
    if(~isfield(Output,'Par'))
        error('Output must have field Par or RecoPar.')
    end
    Output.RecoPar = Output.Par;
end
Output.RecoPar.DataSize = [size_MultiDims(Output.OutTraj.GM,[3 4]) Output.RecoPar.nPartEnc Output.RecoPar.nSLC ...
                           Output.RecoPar.vecSize Output.RecoPar.total_channel_no_measured Output.RecoPar.nRep];

if(~isfield_recursive(Output,'RecoPar.nPartEnc_Meas'))
    Output.RecoPar.nPartEnc_Meas = Output.RecoPar.nPartEnc;
end
                       
if(Output.Par.Hamming_flag && Settings.DensComp_flag)
    fprintf('\nWarning in op_ReadAndReco3DCRTData: Data were Hamming measured, but DensComp was on. Turned DensComp off in Reco.\n')
    Settings.DensComp_flag = false;
end                       
                       
                       
                      
% Output = supp_FixPars(Output);  % To hard-code/hack parameters for special cases, or to make Parameters consistent between different read-in-methods.


tic
fprintf('\n\nReconstructing data\t\t...')




%% Remove first ADC Points

if(isfield(Settings,'RemoveFirstADCPoints'))
    for ii = 1:numel(Output.Data)
        Siz = size(Output.Data{ii}); 
        Tmpp = permute(Output.Data{ii},[1 5 2 3 4 6]);
        Tmpp = reshape(Tmpp,[Siz(1)*Siz(5)/Output.Par.nTempIntsPerAngInt(ii) Output.Par.nTempIntsPerAngInt(ii) Siz(2) Siz(3) Siz(4)]); 
        Start = Settings.RemoveFirstADCPoints{ii};
        RmFIDPts = ceil(Start/Siz(1));
        End = size(Tmpp,1) - (RmFIDPts*Siz(1) - Start + 1);
        Tmpp = Tmpp(Start:End,:,:,:,:,:);
        NewSiz = zeros([1 numel(Siz)]); NewSiz(5) = RmFIDPts*Output.Par.nTempIntsPerAngInt(ii); NewSiz = Siz - NewSiz;
        NewSiz2 = cat(2,NewSiz([1 5 2 3 4]),NewSiz(6:end));
        Tmpp = reshape(Tmpp,NewSiz2);
        Tmpp = permute(Tmpp,[1 3 4 5 2 6]);
        Output.Data{ii} = reshape(Tmpp,NewSiz);
        
        Output.Par.DataSize{ii}(5) = NewSiz(5);
    end
    
    % Noise Data
    if(isfield(Output,'NoiseData'))
        for ii = 1:numel(Output.NoiseData)
            Siz = size(Output.NoiseData{ii}); 
            Tmpp = permute(Output.NoiseData{ii},[1 5 2 3 4 6]);
            Tmpp = reshape(Tmpp,[Siz(1)*Siz(5)/Output.Par.nTempIntsPerAngInt(ii) Output.Par.nTempIntsPerAngInt(ii) Siz(2) Siz(3) Siz(4)]); 
            Start = Settings.RemoveFirstADCPoints{ii};
            RmFIDPts = ceil(Start/Siz(1));
            End = size(Tmpp,1) - (RmFIDPts*Siz(1) - Start + 1);
            Tmpp = Tmpp(Start:End,:,:,:,:,:);
            NewSiz = zeros([1 numel(Siz)]); NewSiz(5) = RmFIDPts*Output.Par.nTempIntsPerAngInt(ii); NewSiz = Siz - NewSiz;
            NewSiz2 = cat(2,NewSiz([1 5 2 3 4]),NewSiz(6:end));
            Tmpp = reshape(Tmpp,NewSiz2);
            Tmpp = permute(Tmpp,[1 3 4 5 2 6]);
            Output.NoiseData{ii} = reshape(Tmpp,NewSiz);
        end    
    end
    
    
    % Find out Min vecSize for all Temporal Interleaves. Cut all FIDs to the same vectorSize
    tmp = cat(1,Output.Par.DataSize{:});
    minvecSize = min(tmp(:,5));
    for ii = 1:numel(Output.Data)
        Output.Data{ii} = Output.Data{ii}(:,:,:,:,1:minvecSize);
        Output.Par.DataSize{ii}(5) = minvecSize;
        
        if(isfield(Output,'NoiseData'))
            Output.NoiseData{ii} = Output.NoiseData{ii}(:,:,:,:,1:minvecSize);
            Output.Par.DataSize{ii}(5) = minvecSize;            
        end
        
    end
    
    Output.RecoPar.vecSize = NewSiz(5);
    Output.Par.vecSize = NewSiz(5);
    Output.RecoPar.DataSize(5) = NewSiz(5);
    clear ii Siz Tmpp Start RmFIDPts End NewSiz;
end




%% FOV SHIFTs Add correct phaseses to the data and shift the FOV to the image center.
LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);

% In Plane Rotation?
Normal2=[1 0 0];
Normal1 = [cos(Output.RecoPar.InPlaneRotation) sin(Output.RecoPar.InPlaneRotation) 0];
v=vrrotvec(Normal1,Normal2);
Rot2=vrrotvec2mat(v);
Rot=Rot2*Rot;


AllFields = fieldnamesr(Output.RecoSteps);
ConjFields = AllFields(~cellfun(@isempty,regexp(AllFields,'ConjInkSpace_flag|ConjIniSpace_flag|ConjFlag')));
ConjNumber = 0; for ii=1:numel(ConjFields); ConjNumber = ConjNumber + eval(['Output.RecoSteps.' ConjFields{ii}]);end
ConjFlag = mod(ConjNumber,2);
if(ConjFlag)
   ConjSign = -1;
else
  ConjSign = 1; 
end


PRS=Rot*LPH';
PRS = PRS + reshape(Settings.ExtraFoVShiftInPatCoordSys,size(PRS));
FOVShift = cellfun( @(x) transpose(exp(ConjSign*1i*x(1,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)), Output.InTraj.GM , 'uni', false);

FOVShift2 = cellfun( @(x) transpose(exp(ConjSign*1i*x(2,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase)),Output.InTraj.GM,'uni',false);

Output.Data = cellfun( @(x,y,z) x.*(y.*z),Output.Data, FOVShift,FOVShift2,'uni',false);


%% Conj in Beginning

if(Settings.ConjInkSpace_flag)
    Output.Data = cellfun(@conj,Output.Data,'uni',0);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData{1}) > 1)
        Output.NoiseData = cellfun(@conj,Output.NoiseData,'uni',0);
    end
end

%% Try to correct for "tilted trajectory" by doing phaseroll

% Save the data for reconstructing the pseudo-pcg case
TiltTrajMat = cell([1 Output.RecoPar.nAngInts]);
if(Settings.Phaseroll_flag)

    nTI = Output.RecoPar.nTempIntsPerAngInt;
    vs = Output.RecoPar.vecSize;
    ns = Output.RecoPar.TrajPts;
    nc = Output.RecoPar.nAngInts;
    nrew = Output.RecoPar.RewPts;
    ncha = size(Output.Data,6);

    
    for ii = 1:nc
%         timeoffset = 0:(ns(ii)-1);
%         timeoffset = repmat(transpose(timeoffset),[1 1 vs]);
%         Freq = ((0:vs-1)/vs-0.5)/(nrew + ns(ii)).*nTI(ii);
%         Freq = myrepmat(Freq,size(timeoffset));    
% 
%         % This comes from:
%         % timeoffset = (0:(ns-1))*Output.RecoPar.ADC_Dt/10^6;
%         % sBW = nTI/((nrew + ns)*Output.RecoPar.ADC_Dt/10^6);
%         % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
%         % Output.RecoPar.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
%         % the rest is basically the same, (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to sBW*((0:vs-1)/vs-0.5) and the other constants are
%         % the same anyway
        
        
        
        % Symmetrized frequencies. Here I assume that my spectrum is symmetric around 0. 
        % E.g. SBW = 20, vs = 20. Can make spectrum go from -10:1:9, or from -9:1:10, or from -9.5:1:9.5. Here I chose the last option. 
        % Before that change, I assumed the first option.
        timeoffset = 0:(ns(ii)-1);
        timeoffset = repmat(transpose(timeoffset),[1 1 vs]);
        Freq = ((0:vs-1)/vs+0.5*(1/vs-1))/(nrew + ns(ii)).*nTI(ii);
        Freq = myrepmat(Freq,size(timeoffset),[0 1 2]);    

        % This comes from:
        % timeoffset = (0:(ns-1))*Output.RecoPar.ADC_Dt/10^6;
        % sBW = nTI/((nrew + ns)*Output.RecoPar.ADC_Dt/10^6);   Why the nTI? Bc if we have 2 TIs it takes us double to finish the circle, but we merge together the FIDs of two TIs and so have the same effective SBW
        % Freq = -sBW/2 + SBW/(2vs) : sBW/vs : (sBW/2 - sBW/(2vs)) = SBW/vs* (-vs/2+0.5 : vs/2-0.5) = SBW/vs* [(0:vs-1) -vs/2+0.5] = 
        % = SBW* [(0:vs-1)/vs - 0.5+0.5/vs] = SBW* [(0:vs-1)/vs + 0.5*(1/vs - 1)] = 
        % Output.RecoPar.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted.
        
        
%         % Lukis pharoll
%         vs = vs*2;
%         PhaseRollDummy = Output.Data{ii};
%         PhaseRollDummy = cat(5,zeros(size(PhaseRollDummy(:,:,:,:,end:-1:1,:))),PhaseRollDummy);
%         timeoffset=nTI(ii)*(((1:ns(ii))-1)/ns(ii)-0); % 3 because of 3 TI. if 2 TI than 2 :)
%         phase(:,:,ii)=-2*1i*pi*timeoffset'*((0:vs-1)/vs-0.5);
%         phasecorr = reshape(exp(phase(:,:,ii)),[ns(ii) 1 1 vs]);
% 
%         kfaxis=squeeze_single_dim(PhaseRollDummy(:,:,:,1,:,:),4);
%         kfaxis=fftshift(ifft(fftshift(kfaxis,4),[],4),4);
%         %kfaxis=circshift(kfaxis,[0 0 -freqshift]);
%         kfaxis=kfaxis.*(conj(phasecorr));
%         %kfaxis=circshift(kfaxis,[0 0 freqshift]);
%         kfaxis=fftshift(fft(fftshift(kfaxis,4),[],4),4);
%         
%         kfaxis=kfaxis(:,:,:,vs/2+1:end,:); 
%         Output.Data{ii}=reshape(kfaxis,size(Output.Data{ii}));
%   
%         TiltTrajMat{ii} = phasecorr;
%         clear kfaxis PhaseRollDummy timeoffset phasecorr phase phasecorr        
%         vs = vs/2;
        
        
        

        phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    % Freq(:,2*end/3)
        Sizzy = size(Output.Data{ii});  Sizzy = cat(2,Sizzy,ones([1 5-numel(Sizzy)]));    % Output.RecoPar.DataSize has wrong size for a short time during the reco. It's alrdy set to the output size
        phasecorr = reshape(phasecorr,[Sizzy(1:2) 1 Sizzy(4:5) 1]); clear Sizzy
    %     phasecorr = myrepmat(phasecorr,size(Output.Data));
    %     phasecorr = conj(phasecorr);
    %     bla = Output.Data(:,:,:,:,1,:,:);
        Output.Data{ii} = fft(fftshift(conj(phasecorr).*fftshift(ifft(Output.Data{ii},[],5),5),5),[],5);
        if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
            Output.NoiseData{ii} = fft(fftshift(conj(phasecorr).*fftshift(ifft(Output.NoiseData{ii},[],5),5),5),[],5);
        end
    %     Output.Data = conj(phasecorr).*Output.Data;     % For correcting only one constant frequency (e.g. at 3ppm which is about the metabo region)

    %     Output.Data(:,:,:,:,1,:,:) = bla;

        TiltTrajMat{ii} = reshape(phasecorr(:,:,1,1,:,1),[Output.RecoPar.TrajPts(ii) 1 Output.RecoPar.vecSize]);  
    end
    
else
    
    for ii = 1:Output.RecoPar.nAngInts
        TiltTrajMat{ii} = ones([Output.RecoPar.TrajPts(ii) Output.RecoPar.vecSize]);
    end
end

if(nargout > 1)
    AdditionalOut.TiltTrajMat = TiltTrajMat;
end


%% Zerofill Circles

% for CurCrcl = 1:Output.RecoPar.nAngInts
% %     NoOfZFPartitions = size(Output.Data{1},3) - size(Output.Data{CurCrcl},3);
%     Output.Data{CurCrcl} = ZerofillOrCutkSpace(Output.Data{CurCrcl},[size_MultiDims(Output.Data{CurCrcl},1:2) size(Output.Data{1},3) size_MultiDims(Output.Data{CurCrcl},4:5)],0);
% end
% Output.RecoPar.nPartEncsPerAngInt(:) = 15;
% Output.RecoPar.AngIntsPerPartEnc = true(size(Output.RecoPar.AngIntsPerPartEnc));


%% Change cell-dependency of nAngInt --> nPart (each Partition one cell-element)
% from {nAngInt}(nTrajPoints x 1 x nPart x nSlc x nTempInt*vecSize x nCha) --> {nPart}(nAngInt*nTrajPoints x Rest)

if(strcmpi(Output.Par.AssumedSequence,'ViennaCRT'))
    Output = op_ChangeCellDependencyFor3DCRT(Output);
elseif(strcmpi(Output.Par.AssumedSequence,'AntoinesEccentricOrRosette'))
    Output = op_ChangeCellDependencyFor3DEccentric(Output);
end


%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape


sft2_Oper = cell([1 numel(Output.InTraj.GM)]);
for ii = 1:numel(Output.InTraj.GM)
    sft2_Oper{ii} = single(sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*Output.RecoPar.DataSize(1)),transpose(Output.InTraj.GM{ii}(:,:)),1));
end


% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(Output.RecoPar.DataSize(1:2)),[1 2],[1 1 1 Output.RecoPar.DataSize(1)/2-1],1); 
    FoVMask = FoVMask(:);
    for ii = 1:numel(Output.InTraj.GM)
        sft2_Oper{ii}(:,~logical(FoVMask(:))) = 0;
    end
    clear FoVMask;
end

%% Calculate B0-Correction of Spiral Data in Spatial Domain
if(Settings.Correct4SpatialB0_flag)
    t = zeros([size(sft2_Oper,1) 1]);
    CurPt = 1;
    for ii = 1:Output.RecoPar.nAngInts
        t(CurPt:CurPt+Output.RecoPar.TrajPts(ii)-1)   = (0:Output.RecoPar.TrajPts(ii)-1)*Output.RecoPar.ADCdtPerAngInt_ns(ii)/10^9;
        CurPt = CurPt + Output.RecoPar.TrajPts(ii);
    end
    CurB0 = imresize(B0.B0Map,Output.RecoPar.DataSize(1:2));    
    if(isfield(B0,'Mask'))
        Mask = imresize(MaskShrinkOrGrow(B0.Mask,2,0,1),Output.RecoPar.DataSize(1:2),'nearest');
        CurB0 = CurB0 .* Mask;
    end
    
    B0CorrMat_Spatial = exp(transpose(-2*pi*1i*CurB0(:)) .* t);

    % SpSpice.Spi2Cart.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;
    % SpSpice.Spi2Cart_NoTiltCorr.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;

    sft2_Oper = sft2_Oper .* B0CorrMat_Spatial;
    
end




%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(Settings.DensComp_flag)
    [Output,Dummy] = op_CalcAndApplyDensComp(Output,sft2_Oper,Settings.DensComp);
    AdditionalOut.DCFPreG = Dummy.DCFPreG; clear Dummy;
end
    
        


    

%% Apply sft2-Operator (Fourier-Transform from spiral k-space --> Cartesian image 
    
Temp = zeros(Output.RecoPar.DataSize,'single');
Temp2 = zeros(Output.RecoPar.DataSize,'single');
for CurPartEnc = 1:Output.RecoPar.nPartEnc_Meas
    
    if(numel(Output.Data) ~= numel(sft2_Oper))      % If for each partition we have identical trajectories, just maybe fewer
        Part = 1;
    else                                          % If each partition has unique trajectories
        Part = CurPartEnc;
    end
    
    Temp(:,:,CurPartEnc,:,:,:,:) = reshape(sft2_Oper{Part}(1:size(Output.Data{CurPartEnc},1),:)' * Output.Data{CurPartEnc} * size(Output.OutTraj.GM(:,:),2),[Output.RecoPar.DataSize(1:2) 1 Output.RecoPar.DataSize(4:end)]);
    Output.Data{CurPartEnc} = [];
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) >= 1) %SB22
            Temp2(:,:,CurPartEnc,:,:,:,:) = reshape(sft2_Oper{Part}(1:size(Output.NoiseData{CurPartEnc},1),:)' * Output.NoiseData{CurPartEnc} * size(Output.OutTraj.GM(:,:),2),[Output.RecoPar.DataSize(1:2) 1 Output.RecoPar.DataSize(4:end)]);
            Output.NoiseData{CurPartEnc} = [];
    end
end
Output.Data = Temp; 
if(isfield(Output,'NoiseData') && numel(Output.NoiseData) >= 1) %SB22
    Output.NoiseData = Temp2;
end
clear Temp Temp2


% Remove fov_overgrid
Output.RecoPar.DataSize(1:2) = Output.RecoPar.DataSize(1:2)/Output.RecoPar.fov_overgrid;
bla = ([size(Output.Data,1) size(Output.Data,2)] - Output.RecoPar.DataSize(1:2))/2+1;
Output.Data = Output.Data(bla(1):bla(1)+Output.RecoPar.DataSize(1)-1,bla(1):bla(2)+Output.RecoPar.DataSize(2)-1,:,:,:,:);
if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
    Output.NoiseData = Output.NoiseData(bla(1):bla(1)+Output.RecoPar.DataSize(1)-1,bla(1):bla(2)+Output.RecoPar.DataSize(2)-1,:,:,:,:);
end

if(nargout > 1)
    AdditionalOut.sft2_Oper = sft2_Oper;
end


%% Perform Reconstruction in Slice and z-dimension

Size = size(Output.Data); Size = cat(2,Size,ones([1 6-numel(Size)]));

if(Size(3) > 1 && Settings.PerformZFFT_flag)
    
    
    
    % z-FT needs to be sFT bc its non-Cartesian
    if(Output.Par.Hamming_flag)
        
        
        kSpaceVector_Normalized = Output.Par.kzPositions/max(abs(Output.Par.kzPositions))/2;
        sft1_Oper = sft1_Operator(kSpaceVector_Normalized,(-floor(Output.Par.nPartEnc/2):1:(ceil(Output.Par.nPartEnc/2)-1)),0); 
        DCF = 1;

        Output.Data = permute(Output.Data,[3 1 2 4 5 6 7]);
        IntermediateSize = size(Output.Data); IntermediateSize(1) = size(sft1_Oper,1);
        Output.Data = reshape(Output.Data,[size(Output.Data,1) numel(Output.Data)/size(Output.Data,1)]);
        Output.Data = reshape(sft1_Oper * (Output.Data.*DCF),IntermediateSize); 
        Output.Data = permute(Output.Data,[2 3 1 4 5 6 7]);
        
        if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
            Output.NoiseData = permute(Output.NoiseData,[3 1 2 4 5 6 7]);
            IntermediateSize = size(Output.NoiseData); IntermediateSize(1) = size(sft1_Oper,1);
            Output.NoiseData = reshape(Output.NoiseData,[size(Output.NoiseData,1) numel(Output.NoiseData)/size(Output.NoiseData,1)]);
            Output.NoiseData = reshape(sft1_Oper * (Output.NoiseData.*DCF),IntermediateSize); 
            Output.NoiseData = permute(Output.NoiseData,[2 3 1 4 5 6 7]);       
        end
       
        
    else
    
        Output.Data = FFTOfMRIData(Output.Data,0,3,0,1,0);
        if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
            Output.NoiseData = FFTOfMRIData(Output.NoiseData,0,3,0,1,0);
        end
    end
    
     
end
% Some special slice-reco necessary? E.g. Hadamard decoding. 
% % if(Size(4) > 1)
% %     
% % end

% % Merge slice and z-dimensions. Normally they are exclusive anyway...
% Output.Data = reshape(Output.Data, [Size(1:2) prod(Size(3:4)) Size(5:end)]); 
% if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
%     Output.NoiseData = reshape(Output.NoiseData, [Size(1:2) prod(Size(3:4)) Size(5:end)]); 
% end


%% Conj at End

if(Settings.ConjIniSpace_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Output.NoiseData = conj(Output.NoiseData);
    end
end


%% Flip left right

if(Settings.FlipDim_flag)
    Output.Data = flip(Output.Data,Settings.FlipDim);
%     Output.Data = circshift(Output.Data,[1 -1 0 0 0 0 0]);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Output.NoiseData = flip(Output.NoiseData,Settings.FlipDim);
%         Output.NoiseData = circshift(Output.NoiseData,[1 -1 0 0 0 0 0]);        
    end
end


%% Postparations

Output = supp_UpdateRecoSteps(Output,Settings);

fprintf('\n\t\t\t\t...\ttook\t%10.6f seconds',toc)


