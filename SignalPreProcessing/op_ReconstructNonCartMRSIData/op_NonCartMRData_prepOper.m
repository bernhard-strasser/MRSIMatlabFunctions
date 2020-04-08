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
Operators.InDataSize = size(Output.Data);

Operators.OutDataSize = [size_MultiDims(Output.OutTraj.GM,[3 4]) Output.RecoPar.nPartEnc*Output.RecoPar.nSLC Output.RecoPar.vecSize];
if(isfield(Output.Par,'TimeUndersamplFactor'))
    Operators.OutDataSize(end) = Operators.OutDataSize(end)*Output.Par.TimeUndersamplFactor;
end
Output.RecoPar.DataSize = Operators.OutDataSize;
% Output.Par.total_channel_no_measured: Can we somehow find out if we will do a coil combination in our reco or not?


%% FoV-Shift Operator

LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);
PRS=Rot*LPH';
Operators.FoVShift = squeeze(exp(1i*Output.InTraj.GM(1,:,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)); 
Operators.FoVShift = Operators.FoVShift .* squeeze(exp(1i*Output.InTraj.GM(2,:,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase));
clear LPH Normal1 Normal2 v Rot PRS 



%% Sampling Operator

if(~isfield(AdditionalIn,'SamplingOperator'))
    Operators.SamplingOperator = 1;
else
    Operators.SamplingOperator = AdditionalIn.SamplingOperator;
    SizeSamp = size(Operators.SamplingOperator);
    UnequalSizesInd = find(Operators.InDataSize(1:end-1) ~= SizeSamp,1);
    if(~isempty(UnequalSizesInd))
        st = dbstack;
        FunName = st(1).name;
        fprintf('\nWarning in %s: Size of given SamplingOperator (size: %s)\ndoes not match InDataSize (%s). Cut SamplingOperator.',FunName,sprintf('%d ',SizeSamp),sprintf('%d ',Operators.InDataSize))
    end
    Operators.SamplingOperator = Zerofilling_Spectral(Operators.SamplingOperator,Operators.InDataSize(1:end-1),0);
end


%% "Phaseroll" Operator
% (The different k-space points were not acquired at same time points, but at time points (t0,t1,t2,..., tN). Therefore, different chemical species with off-resonance
% frequencies (f1,f2,...) will acquire different linear phases along the trajectory, e.g. f1*(t0,t1,t2,...,tN) or f2*(t0,t1,t2,...,tN). This linear phase affects the images,
% for spirals this would be a kind of blurring, which depends on the off-resonance frequencies. For EPSI, it would be shifts, for circles it would be rotations. 
% Try to undo this effect by Fourier transforming each k-space point to spectral domain, multiplying the inverse of this phase for each k-space point (depending on
% frequencies and k-space points (=time points), and Fourier transforming back to time domain.)

if(Settings.Phaseroll_flag)

    nTI = Output.RecoPar.nTempInt;
    vs = Output.RecoPar.vecSize;
    ns = Output.RecoPar.TrajPts;
    nc = Output.RecoPar.nAngInts;
    nrew = Output.RecoPar.RewPts;
    ncha = Operators.InDataSize(6);

    timeoffset = 0:(ns-1);
    timeoffset = repmat(transpose(timeoffset),[1 vs]);
    Freq = ((0:vs-1)/vs-0.5)/(nrew + ns)*nTI;
    Freq = repmat(Freq,[ns 1]);    

    % This comes from:
    % timeoffset = (0:(ns-1))*Output.Par.ADC_Dt/10^6;
    % sBW = nTI/((nrew + ns)*Output.Par.ADC_Dt/10^6);
    % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
    % Output.Par.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
    % the rest is basically the same (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to ((0:vs-1)/vs-0.5), and the other constants are
    % the same anyway
    
    Operators.TiltTrajMat = exp(-2*1i*pi*timeoffset .* Freq);    
    Operators.TiltTrajMat = reshape(Operators.TiltTrajMat,[Operators.InDataSize(1) 1 1 1 Operators.InDataSize(5)]); % [TrajPts nAI nPart nSlc vecSize]
       
else
    
	Operators.TiltTrajMat = 1;

end



%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])

Operators.sft2_Oper = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*Operators.OutDataSize(1)),transpose(Output.InTraj.GM(:,:)),1);

% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(Operators.OutDataSize(1:2)),[1 2],[1 1 1 Operators.OutDataSize(1)/2-1],1); 
    FoVMask = FoVMask(:);
    Operators.sft2_Oper(:,~logical(FoVMask(:))) = 0;
    clear FoVMask;
end


%% Calculate B0-Correction of Spiral Data in Spatial Domain
if(Settings.Correct4SpatialB0_flag)
    t   = (0:Output.RecoPar.TrajPts-1)*Output.RecoPar.ADC_dt/10^9;
    t = repmat(t,[1 1 Output.RecoPar.nAngInts]); t = t(:);
    CurB0 = imresize(AdditionalIn.B0.B0Map,Operators.OutDataSize(1:2));    
    if(isfield(AdditionalIn.B0,'Mask'))
        Mask = imresize(MaskShrinkOrGrow(AdditionalIn.B0.Mask,2,0,1),Operators.OutDataSize(1:2),'nearest');
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

    time   = (0:Operators.InDataSize(5)-1)*Output.RecoPar.Dwelltimes(1)/10^9; time = reshape(time,[1 1 1 numel(time)]);
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
    Operators.DCFPreG = 1;
end


%% Calculate SensitivityMap

Settss.ResizeMethod = 'Imresize';
Settss.ScalingMethod = 'UniformSensitivity';
Dummy2.Data = ones([Operators.OutDataSize Output.Par.total_channel_no_measured]); 
Dummy2.RecoPar = Output.RecoPar; Dummy2.RecoPar.DataSize = size(Dummy2.Data);  %Dummy2.Par = Output.Par; 
[DeleteMe, AddOut] = op_CoilCombineData(Dummy2,AdditionalIn.SensMap,Settss);
Operators.SensMap = AddOut.CoilWeightMap(:,:,:,1,:) .* AddOut.Scaling;
clear DeleteMe AddOut Settss;

Operators.SensMap(isinf(Operators.SensMap) | isnan(Operators.SensMap)) = 0;



