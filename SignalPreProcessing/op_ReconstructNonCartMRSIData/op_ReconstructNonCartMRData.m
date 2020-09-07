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
if(~isfield(Settings,'FlipDim_flag'))
    Settings.FlipDim_flag = false;
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
Output.RecoPar.DataSize = [size_MultiDims(Output.OutTraj.GM,[3 4]) Output.RecoPar.nPartEnc Output.RecoPar.nSLC ...
                           Output.RecoPar.vecSize Output.RecoPar.total_channel_no_measured];


                      
% Output = supp_FixPars(Output);  % To hard-code/hack parameters for special cases, or to make Parameters consistent between different read-in-methods.


tic
fprintf('\n\nReconstructing data\t\t...')


%% FOV SHIFTs Add correct phaseses to the data and shift the FOV to the image center.
LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);

% In Plane Rotation?
% Normal2=[1 0 0];
% Normal1 = [cos(InPlaneRot) sin(InPlaneRot) 0];
% v=vrrotvec(Normal1,Normal2);
% Rot2=vrrotvec2mat(v);
% Rot=Rot2*Rot;

PRS=Rot*LPH';
FOVShift = cellfun( @(x) transpose(exp(1i*x(1,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)), Output.InTraj.GM , 'uni', false);

FOVShift2 = cellfun( @(x) transpose(exp(1i*x(2,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase)),Output.InTraj.GM,'uni',false);

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



%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nTrajPoints x nAngInt x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])
Output.Data = cat(1,Output.Data{:});
SizeData_k = size(Output.Data); SizeData_k = cat(2,SizeData_k,ones([1 5-numel(SizeData_k)]));
Output.Data = Output.Data(:,:);
if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
    Output.NoiseData = cat(1,Output.NoiseData{:});
    Output.NoiseData = Output.NoiseData(:,:);   
end


sft2_Oper = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*Output.RecoPar.DataSize(1)),transpose(cat(2,Output.InTraj.GM{:})),1);

% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(Output.RecoPar.DataSize(1:2)),[1 2],[1 1 1 Output.RecoPar.DataSize(1)/2-1],1); 
    FoVMask = FoVMask(:);
    sft2_Oper(:,~logical(FoVMask(:))) = 0;
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
    
Output.Data = sft2_Oper' * Output.Data * size(Output.OutTraj.GM(:,:),2);
Output.Data = reshape(Output.Data,[Output.RecoPar.DataSize(1:2) SizeData_k(3:end)]);
if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
    Output.NoiseData = sft2_Oper' * Output.NoiseData * size(Output.OutTraj.GM(:,:),2);
    Output.NoiseData = reshape(Output.NoiseData,[Output.RecoPar.DataSize(1:2) SizeData_k(3:end)]);
end

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

Size = size(Output.Data);

if(Size(3) > 1 && Settings.PerformZFFT_flag)
    Output.Data = FFTOfMRIData(Output.Data,0,3,0,1,0);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Output.NoiseData = FFTOfMRIData(Output.NoiseData,0,3,0,1,0);
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
    Output.Data = flip(Output.Data,1);
%     Output.Data = circshift(Output.Data,[1 -1 0 0 0 0 0]);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Output.NoiseData = flip(Output.NoiseData,1);
%         Output.NoiseData = circshift(Output.NoiseData,[1 -1 0 0 0 0 0]);        
    end
end


%% Postparations

Output = supp_UpdateRecoSteps(Output,Settings);

fprintf('\n\t\t\t\t...\ttook\t%10.6f seconds',toc)


