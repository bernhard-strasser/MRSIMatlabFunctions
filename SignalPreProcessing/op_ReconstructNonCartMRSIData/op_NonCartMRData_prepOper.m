function Operators = op_NonCartMRData_prepOper(Output,AdditionalIn,Settings)
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
% [Output, AdditionalOut] = op_NonCartMRData_prepOper(Output,AdditionalIn,Settings)
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
if(~isfield_recursive(Settings,'DensComp.AutoScale_flag'))
    Settings.DensComp.AutoScale_flag = false;    
end
if(~isfield(Settings.DensComp,'Normalize_flag'))
    Settings.DensComp.Normalize_flag = false;    
end
if(~isfield(Settings,'Correct4SpatialB0_flag'))
   Settings.Correct4SpatialB0_flag = false;    
end
if(~isfield(Settings,'Correct4SpectralB0_flag'))
   Settings.Correct4SpectralB0_flag = false;    
end
if(~isfield(Settings,'CircularSFTFoV_flag'))
   Settings.CircularSFTFoV_flag = false;    
end

if(~isfield(Output,'RecoPar'))
    if(~isfield(Output,'Par'))
        error('Output must have field Par or RecoPar.')
    end
    Output.RecoPar = Output.Par;
end
Operators.InDataSize = Output.RecoPar.DataSize;

Operators.OutDataSize = [size_MultiDims(Output.OutTraj.GM,[3 4]) Output.RecoPar.nPartEnc*Output.RecoPar.nSLC Output.RecoPar.vecSize];
% if(isfield(Output.Par,'TimeUndersamplFactor'))
%     Operators.OutDataSize(end) = Operators.OutDataSize(end)*Output.Par.TimeUndersamplFactor;
% end
Output.RecoPar.DataSize = Operators.OutDataSize;
% Output.Par.total_channel_no_measured: Can we somehow find out if we will do a coil combination in our reco or not?

Output = supp_FixPars(Output);  % To hard-code/hack parameters for special cases, or to make Parameters consistent between different read-in-methods.



%% FoV-Shift Operator

% LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
% Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
% Normal2=[0 0 1];
% v=vrrotvec(Normal1,Normal2);
% Rot=vrrotvec2mat(v);
% PRS=Rot*LPH';
% Operators.FoVShift = squeeze(exp(1i*Output.InTraj.GM(1,:,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)); 
% Operators.FoVShift = Operators.FoVShift .* squeeze(exp(1i*Output.InTraj.GM(2,:,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase));
% clear LPH Normal1 Normal2 v Rot PRS 


LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);
PRS=Rot*LPH';
FOVShift = cellfun( @(x) transpose(exp(1i*x(1,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)), Output.InTraj.GM , 'uni', false);
FOVShift2 = cellfun( @(x) transpose(exp(1i*x(2,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase)),Output.InTraj.GM,'uni',false);
Operators.FoVShift = cellfun(@times,FOVShift,FOVShift2,'uni',0);




%% Sampling Operator

if(~isfield(AdditionalIn,'SamplingOperator'))
    Operators.SamplingOperator = 1;
else
    Operators.SamplingOperator = AdditionalIn.SamplingOperator;
    SizeSamp = cellfun(@size,Operators.SamplingOperator,'uni',false);
    
    test = cat(1,Operators.InDataSize{:}); 
    
    [UnequalSizesInd_AngIntDim,UnequalSizesInd_VecDim] = find(test(:,1:5) ~= cat(1,SizeSamp{:}));
    if(~isempty(UnequalSizesInd_AngIntDim))
        st = dbstack;
        FunName = st(1).name;
        fprintf('\nWarning in %s: Size of given SamplingOperator (size: %s)\ndoes not match InDataSize (%s). Cut SamplingOperator.',...
        FunName,sprintf(['[' repmat('%d ',[1 numel(SizeSamp{UnequalSizesInd_AngIntDim(1)})]) '], '],(cat(3,SizeSamp{UnequalSizesInd_AngIntDim}))),...
        sprintf(['[' repmat('%d ',[1 numel(Operators.InDataSize{UnequalSizesInd_AngIntDim(1)})]) '], '],(cat(3,Operators.InDataSize{UnequalSizesInd_AngIntDim}))))

    end
    for CurAI = UnequalSizesInd_AngIntDim'
        Operators.SamplingOperator{CurAI} = Zerofilling_Spectral(Operators.SamplingOperator{CurAI},Operators.InDataSize{CurAI}(1:end-1),0);
    end
end


%% "Phaseroll" Operator
% (The different k-space points were not acquired at same time points, but at time points (t0,t1,t2,..., tN). Therefore, different chemical species with off-resonance
% frequencies (f1,f2,...) will acquire different linear phases along the trajectory, e.g. f1*(t0,t1,t2,...,tN) or f2*(t0,t1,t2,...,tN). This linear phase affects the images,
% for spirals this would be a kind of blurring, which depends on the off-resonance frequencies. For EPSI, it would be shifts, for circles it would be rotations. 
% Try to undo this effect by Fourier transforming each k-space point to spectral domain, multiplying the inverse of this phase for each k-space point (depending on
% frequencies and k-space points (=time points), and Fourier transforming back to time domain.)

Operators.TiltTrajMat = cell([1 Output.RecoPar.nAngInts]);
if(Settings.Phaseroll_flag)

    nTI = Output.RecoPar.nTempIntsPerAngInt;
    vs = Output.RecoPar.vecSize;
    ns = Output.RecoPar.TrajPts;
    nc = Output.RecoPar.nAngInts;
    nrew = Output.RecoPar.RewPts;
    ncha = size(Output.Data,6);

    
    for ii = 1:nc
        timeoffset = 0:(ns(ii)-1);
        timeoffset = repmat(transpose(timeoffset),[1 1 vs]);
        Freq = ((0:vs-1)/vs-0.5)/(nrew + ns(ii)).*nTI(ii);
        Freq = myrepmat(Freq,size(timeoffset));    

        % This comes from:
        % timeoffset = (0:(ns-1))*Output.RecoPar.ADC_Dt/10^6;
        % sBW = nTI/((nrew + ns)*Output.RecoPar.ADC_Dt/10^6);
        % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
        % Output.RecoPar.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
        % the rest is basically the same (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to ((0:vs-1)/vs-0.5), and the other constants are
        % the same anyway

        phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    % Freq(:,2*end/3)
        Sizzy = size(Output.Data{ii});      % Output.RecoPar.DataSize has wrong size for a short time during the reco. It's alrdy set to the output size
        phasecorr = reshape(phasecorr,[Sizzy(1:5) 1]); clear Sizzy
    %     phasecorr = myrepmat(phasecorr,size(Output.Data));
    %     phasecorr = conj(phasecorr);
    %     bla = Output.Data(:,:,:,:,1,:,:);
        Output.Data{ii} = fft(fftshift(conj(phasecorr).*fftshift(ifft(Output.Data{ii},[],5),5),5),[],5);
        if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
            Output.NoiseData{ii} = fft(fftshift(conj(phasecorr).*fftshift(ifft(Output.NoiseData{ii},[],5),5),5),[],5);
        end
    %     Output.Data = conj(phasecorr).*Output.Data;     % For correcting only one constant frequency (e.g. at 3ppm which is about the metabo region)

    %     Output.Data(:,:,:,:,1,:,:) = bla;

        Operators.TiltTrajMat{ii} = phasecorr(:,:,1,1,:,1);  
    end
    
else
    
    for ii = 1:Output.RecoPar.nAngInts
        Operators.TiltTrajMat{ii} = ones([Output.RecoPar.TrajPts(ii) Output.RecoPar.vecSize]);
    end
end



%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])

Operators.sft2_Oper = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*Operators.OutDataSize(1)),transpose(cat(2,Output.InTraj.GM{:})),1);

% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(Operators.OutDataSize(1:2)),[1 2],[1 1 1 Operators.OutDataSize(1)/2-1],1); 
    FoVMask = FoVMask(:);
    Operators.sft2_Oper(:,~logical(FoVMask(:))) = 0;
    clear FoVMask;
end


%% Calculate B0-Correction of Spiral Data in Spatial Domain
if(Settings.Correct4SpatialB0_flag)
    t = zeros([size(Operators.sft2_Oper,1) 1]);
    CurPt = 1;
    for ii = 1:Output.RecoPar.nAngInts
        t(CurPt:CurPt+Output.RecoPar.TrajPts(ii)-1)   = (0:Output.RecoPar.TrajPts(ii)-1)*Output.RecoPar.ADCdtPerAngInt_ns(ii)/10^9;
        CurPt = CurPt + Output.RecoPar.TrajPts(ii);
    end
    CurB0 = imresize(AdditionalIn.B0.B0Map,Output.RecoPar.DataSize(1:2));    
    if(isfield(AdditionalIn.B0,'Mask'))
        Mask = imresize(MaskShrinkOrGrow(AdditionalIn.B0.Mask,2,0,1),Output.RecoPar.DataSize(1:2),'nearest');
        CurB0 = CurB0 .* Mask;
    end
    
    Operators.B0CorrMat_Spatial = exp(transpose(-2*pi*1i*CurB0(:)) .* t);

    Operators.sft2_Oper = Operators.sft2_Oper .* Operators.B0CorrMat_Spatial;
end


%% Mask

if(isfield_recursive(AdditionalIn,'B0.BrainMask'))
    Mask = imresize(AdditionalIn.B0.BrainMask,Operators.OutDataSize(1:2),'nearest');
    Operators.Mask = Mask;
else
    Mask = ones(Operators.OutDataSize(1:2)); 
end


%% Calculate B0CorrMat_Spec

if(Settings.Correct4SpectralB0_flag)
    % Go to resolution of D1
    CurB0 = imresize(AdditionalIn.B0.B0Map,Operators.OutDataSize(1:2)) ;
    if(exist('Mask','var'))
        CurB0 = CurB0 .* Mask;
    end

    % Round to shift only integer number of points
    HzPerPt = 10^9/Output.RecoPar.Dwelltimes(1) / Output.RecoPar.vecSize;
    % if(Settings.RoundB0ToIntVecPts)
    %     CurB0 = round(CurB0/HzPerPt)*HzPerPt;
    % end

    time   = (0:Operators.OutDataSize(4)-1)*Output.RecoPar.Dwelltimes(1)/10^9; time = reshape(time,[1 1 1 numel(time)]);
    B0CorrMat_Spec = exp(2*pi*1i*CurB0 .* time);


    Operators.B0CorrMat_Spec = B0CorrMat_Spec;
else
    Operators.B0CorrMat_Spec = 1;
end


%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(Settings.DensComp_flag)
    if(isfield(Output,'Data'))
        Output = rmfield(Output,'Data'); 
    end
    if(isfield(Output,'NoiseData'))
        Output = rmfield(Output,'NoiseData');
    end 
    [~,Dummy2] = op_CalcAndApplyDensComp(Output,Operators.sft2_Oper,Settings.DensComp);
    Operators.DCFPreG = Dummy2.DCFPreG;
    clear Dummy;
else
%     Operators.DCFPreG = 1;
    ImOut = Operators.sft2_Oper' * (Operators.sft2_Oper * ones([size(Operators.sft2_Oper,2) 1]) * size(Output.OutTraj.GM(:,:),2));
    ImOut(ImOut == 0) = NaN; 
    Operators.DCFPreG = 1/nanmean(abs(ImOut));
    
end


%% Calculate SensitivityMap

if(isfield(AdditionalIn,'SensMap'))
    Settss.ResizeMethod = 'Imresize';
    Settss.ScalingMethod = 'UniformSensitivity';
    Dummy2.Data = ones([Operators.OutDataSize Output.Par.total_channel_no_measured]); 
    Dummy2.RecoPar = Output.RecoPar; Dummy2.RecoPar.DataSize = size(Dummy2.Data);  %Dummy2.Par = Output.Par; 
    [DeleteMe, AddOut] = op_CoilCombineData(Dummy2,AdditionalIn.SensMap,Settss);
    Operators.SensMap = AddOut.CoilWeightMap(:,:,:,1,:) .* AddOut.Scaling;
    clear DeleteMe AddOut Settss;

    Operators.SensMap(isinf(Operators.SensMap) | isnan(Operators.SensMap)) = 0;
else
    Operators.SensMap = 1;
end




%%
if(isfield(Settings,'Uncell_flag') && Settings.Uncell_flag)
    dummy = sum(cat(1,Operators.InDataSize{:}));
    Operators.InDataSize = [dummy(1) Operators.InDataSize{1}(2:end)];
    Operators.FoVShift = cat(1,Operators.FoVShift{:});
    Operators.SamplingOperator = cat(1,Operators.SamplingOperator{:});
    Operators.TiltTrajMat = cat(1,Operators.TiltTrajMat{:});
    
end


