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


FunExists = which('fieldnamesr');
if(isempty(FunExists))
    error('Function depends on function ''fieldnamesr''.\nPlease download e.g. from here: https://de.mathworks.com/matlabcentral/fileexchange/33262-get-structure-field-names-in-recursive-manner')
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
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Output.NoiseData = conj(Output.NoiseData);
    end
end

%% Perform Construction in Slice and z-dimension

Size = size(Output.Data);
if(Size(3) > 1 && Settings.PerformZFFT_flag)
    Output.Data = FFTOfMRIData(Output.Data,0,3,1,1,0);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Output.NoiseData = FFTOfMRIData(Output.NoiseData,0,3,1,1,0);
    end
end

%% Calculate & Apply sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])
SizeData_i = size(Output.Data); SizeData_i = cat(2,SizeData_i,ones([1 5-numel(SizeData_i)]));
sft2_Oper = sft2_Operator(transpose(squeeze(Output.OutTraj.GM(:,:))*SizeData_i(1)),transpose(cat(2,Output.InTraj.GM{:})),1);



% Restrict to circular FoV
if(Settings.CircularSFTFoV_flag)
    FoVMask = EllipticalFilter(ones(size_MultiDims(Output.OutTraj.GM,[1 2])),[1 2],[1 1 1 size(Output.OutTraj.GM,1)/2-1],1); 
    FoVMask = FoVMask(:);
    sft2_Oper(:,~logical(FoVMask(:))) = 0;
    clear FoVMask;
end

Output.Data = reshape(Output.Data,[prod(SizeData_i(1:2)) prod(SizeData_i(3:end))]);   
Output.Data = sft2_Oper * Output.Data;  %/ size(Output.OutTraj.GM(:,:),2)

if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
    Output.NoiseData = reshape(Output.NoiseData,[prod(SizeData_i(1:2)) prod(SizeData_i(3:end))]);   
    Output.NoiseData = sft2_Oper * Output.NoiseData;  %/ size(Output.OutTraj.GM(:,:),2)
    Output.NoiseData = reshape(Output.NoiseData,[Output.Par.DataSize(1:2) SizeData_i(3:end)]);
end

if(nargout > 1)
    AdditionalOut.sft2_Oper = sft2_Oper;
end



%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(Settings.DensComp_flag)
    [Output,Dummy] = op_CalcAndApplyDensComp(Output,sft2_Oper,Settings.DensComp);
    AdditionalOut.DCFPreG = Dummy.DCFPreG; clear Dummy;
end


%% Make Cell
CurPt = 1;
for CurAI = 1:Output.RecoPar.nAngInts
    Dummy{CurAI} = reshape(Output.Data(CurPt:CurPt+Output.Par.DataSize{CurAI}(1)-1,:),[Output.Par.DataSize{CurAI}(1:2) SizeData_i(3:end)]);
    CurPt = CurPt+Output.Par.DataSize{CurAI}(1);
end
Output.Data = Dummy; clear Dummy;

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
        Output.Data{ii} = fft(fftshift((phasecorr).*fftshift(ifft(Output.Data{ii},[],5),5),5),[],5);
        if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
            Output.NoiseData{ii} = fft(fftshift((phasecorr).*fftshift(ifft(Output.NoiseData{ii},[],5),5),5),[],5);
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
    Output.Data = cellfun(@conj,Output.Data,'uni',0);
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData{1}) > 1)
        Output.NoiseData = cellfun(@conj,Output.NoiseData,'uni',0);
    end
end


%% FOV SHIFTs Add correct phaseses to the data and shift the FOV to the image center.
LPH=[Output.RecoPar.Pos_Cor Output.RecoPar.Pos_Sag Output.RecoPar.Pos_Tra];
Normal1=[Output.RecoPar.SliceNormalVector_y Output.RecoPar.SliceNormalVector_x Output.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);


AllFields = fieldnamesr(Output.RecoSteps);
ConjFields = AllFields(~cellfun(@isempty,regexp(AllFields,'ConjInkSpace_flag|ConjIniSpace_flag|ConjFlag')));
ConjNumber = 0; for ii=1:numel(ConjFields); ConjNumber = ConjNumber + eval(['Output.RecoSteps.' ConjFields{ii}]);end
ConjFlag = mod(ConjNumber,2);
if(ConjFlag)
   ConjSign = -1;
else
  ConjSign = 1; 
end



% In Plane Rotation?
% Normal2=[1 0 0];
% Normal1 = [cos(InPlaneRot) sin(InPlaneRot) 0];
% v=vrrotvec(Normal1,Normal2);
% Rot2=vrrotvec2mat(v);
% Rot=Rot2*Rot;

PRS=Rot*LPH';
FOVShift = cellfun( @(x) transpose(exp(-ConjSign*1i*x(1,:)/0.5*Output.RecoPar.DataSize(2)*pi*-PRS(2)/Output.RecoPar.FoV_Read)), Output.InTraj.GM , 'uni', false);

FOVShift2 = cellfun( @(x) transpose(exp(-ConjSign*1i*x(2,:)/0.5*Output.RecoPar.DataSize(1)*pi*PRS(1)/Output.RecoPar.FoV_Phase)),Output.InTraj.GM,'uni',false);

Output.Data = cellfun( @(x,y,z) x.*(y.*z),Output.Data, FOVShift,FOVShift2,'uni',false);


%% Postparations

Output.RecoPar.DataSize = cellfun(@size,Output.Data,'uni',0);

Output = supp_UpdateRecoSteps(Output,Settings);



