function [Output, AdditionalOut] = op_ConstructNonCartMRData(Output,AdditionalIn,Settings)
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
% [Output, AdditionalOut] = op_ConstructNonCartMRData(Output,AdditionalIn,Settings)
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


if(~exist('Settings','var'))
   Settings = struct();  
end
if(~isfield(Settings,'Phaseroll_flag'))
   Settings.Phaseroll_flag = true;    
end
if(~isfield(Settings,'DensComp_flag'))
   Settings.DensComp_flag = true;    
end
if(~isfield(Settings,'DensCompAutoScale_flag'))
   Settings.DensCompAutoScale_flag = false;    
end
if(~isfield(Settings,'ConjIniSpace_flag'))
   Settings.ConjIniSpace_flag = false;    
end
if(~isfield(Settings,'ConjInkSpace_flag'))
   Settings.ConjInkSpace_flag = true;    
end
if(~isfield(Settings,'Correct4SpatialB0_flag'))
   Settings.Correct4SpatialB0_flag = false;    
end
if(~isfield(Settings,'CircularSFTFoV_flag'))
   Settings.CircularSFTFoV_flag = false;    
end
if(exist('AdditionalIn','var') && isfield('AdditionalIn','B0'))
    B0 = AdditionalIn.B0;
end
if(~isfield(Output,'RecoPar'))
    if(~isfield(Output,'Par'))
        error('Output must have field Par or RecoPar.')
    end
    Output.RecoPar = Output.Par;
end



%% Conj in Beginning

if(Settings.ConjIniSpace_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end

%% Perform Construction in Slice and z-dimension

Size = size(Output.Data);
if(Size(3) > 1 && Settings.PerformZFFT_flag)
    Output.Data = FFTOfMRIData(Output.Data,0,3,1,1,0);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = FFTOfMRIData(Output.NoiseData,0,3,1,1,0);
    end
end

%% Calculate & Apply sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])
SizeData_i = size(Output.Data); SizeData_i = cat(2,SizeData_i,ones([1 5-numel(SizeData_i)]));
sft2_Op = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*SizeData_i(1)),transpose(Output.InTraj.GM(:,:)),1);

% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(size_MultiDims(Output.OutTraj.GM,[1 2])),[1 2],[1 1 1 size(Output.OutTraj.GM,1)/2-1],1); 
    FoVMask = FoVMask(:);
    sft2_Oper(:,~logical(FoVMask(:))) = 0;
    clear FoVMask;
end

Output.Data = reshape(Output.Data,[prod(SizeData_i(1:2)) prod(SizeData_i(3:end))]);   
Output.Data = sft2_Op * Output.Data;  %/ size(Output.OutTraj.GM(:,:),2)
Output.Data = reshape(Output.Data,[Output.Par.DataSize(1:2) SizeData_i(3:end)]);

if(isfield(Output,'NoiseData'))
    Output.NoiseData = reshape(Output.NoiseData,[prod(SizeData_i(1:2)) prod(SizeData_i(3:end))]);   
    Output.NoiseData = sft2_Op * Output.NoiseData;  %/ size(Output.OutTraj.GM(:,:),2)
    Output.NoiseData = reshape(Output.NoiseData,[Output.Par.DataSize(1:2) SizeData_i(3:end)]);
end

if(nargout > 1)
    AdditionalOut.sft2_Op = sft2_Op;
end



%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(Settings.DensComp_flag)
    [Output,Dummy] = op_CalcAndApplyDensComp(Output,sft2_Oper,Settings.DensComp);
    AdditionalOut.DCFPreG = Dummy.DCFPreG; clear Dummy;
end


%% Try to correct for "tilted trajectory" by doing phaseroll

% Save the data for reconstructing the pseudo-pcg case
if(Settings.Phaseroll_flag)

    nTI = Output.Par.nTempInt;
    vs = Output.Par.vecSize;
    ns = Output.Par.TrajPts;
    nc = Output.Par.nAngInts;
    nrew = Output.Par.RewPts;
    ncha = size(Output.Data,6);

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
    
    phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    
    phasecorr = myrepmat(phasecorr,size(Output.Data));
    
    bla = Output.Data(:,:,:,:,1,:,:);
    Output.Data = fft(fftshift((phasecorr).*fftshift(ifft(Output.Data,[],5),5),5),[],5);
    Output.Data(:,:,:,:,1,:,:) = bla;
    
    if(isfield(Output,'NoiseData'))
        bla = Output.NoiseData(:,:,:,:,1,:,:);        
        Output.NoiseData = fft(fftshift((phasecorr).*fftshift(ifft(Output.NoiseData,[],5),5),5),[],5);
        Output.NoiseData(:,:,:,:,1,:,:) = bla;
    end    

    TiltTrajMat = phasecorr;   
    
else
    
	TiltTrajMat = ones([Output.Par.TrajPts*Output.Par.nAngInts Output.Par.vecSize]);

end

if(nargout > 1)
    AdditionalOut.TiltTrajMat = TiltTrajMat;
end



%% FOV SHIFTs

% NOT YET IMPLEMENTED

% THATS HOW LUKAS HINGERL DOES IT (CODE FROM LUKAS HINGERL):
% %Inplane FOV shift
% kspace_max=sqrt(kx(1,end).^2+ky(1,end).^2);
% yIPFOVshift=-ReadInInfo.General.Ascconv.PosVOI_Sag; %this should be right!
% xIPFOVshift=-ReadInInfo.General.Ascconv.PosVOI_Cor;
% zIPFOVshift=ReadInInfo.General.Ascconv.PosVOI_Tra;
% FOV=ReadInInfo.General.Ascconv.FoV_Phase;
% FOVz=ReadInInfo.General.Ascconv.FoV_Partition;
% 
% 
% for circles=1:nc
%     for samplepoint=1:ns
%         csi_k(:,circles,samplepoint,1,:)=csi_k(:,circles,samplepoint,1,:)*exp(1i*kx(samplepoint,circles)/kspace_max*((REShalbe-0.5)/FOV)*2*pi*xIPFOVshift)*exp(1i*ky(samplepoint,circles)/kspace_max*((REShalbe-0.5)/FOV)*2*pi*yIPFOVshift);
%     end
% end


%% Conj at End

if(Settings.ConjInkSpace_flag)
    Output.Data = conj(Output.Data);
    
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end    
end


%% FOV SHIFTs Add correct phaseses to the data and shift the FOV to the image center.
LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);
PRS=Rot*LPH';
FOVShift = squeeze(exp(-1i*Output.InTraj.GM(1,:,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)); 
FOVShift = FOVShift .* squeeze(exp(-1i*Output.InTraj.GM(2,:,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase));
Output.Data = Output.Data.*FOVShift;


%% Postparations

Output.RecoPar.DataSize = size(Output.Data);

Output = supp_UpdateRecoSteps(Output,Settings);



