function [Output, AdditionalOut] = op_IterReconstructNonCartMRData(Output,AdditionalIn,ModelFunction,Settings)
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
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
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
if(~isfield_recursive(Settings,'DensComp.AutoScale_flag'))
   Settings.DensComp.AutoScale_flag = false;    
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
Output.RecoPar.DataSize = [size_MultiDims(Output.OutTraj.GM,[3 4]) Output.RecoPar.nPartEnc*Output.RecoPar.nSLC Output.RecoPar.vecSize];
% Output.Par.total_channel_no_measured: Can we somehow find out if we will do a coil combination in our reco or not?

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



%% Conj in Beginning

if(Settings.ConjInkSpace_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end


%% Iterative Reconstruction of CSI Data

% Prepare operators
Opers = op_NonCartMRData_prepOper(Output,AdditionalIn,Settings);

AO = @(x) ModelFunction('NoTransj',x,Opers);
AOT = @(x) ModelFunction('Transj',x,Opers);

% [m,n,p] = size(Opers.B0CorrMat_Spec);
maxiter = 20;
lambda = 0.001;%1e-16;%0.001;
eta =20e-2;

% Output.Data = lr_method(Output.Data,AO,AOT,maxiter,lambda,eta,m,n,Settings, Output.RecoPar);
Output.Data = LS_method(Output.Data,AO,AOT,maxiter,eta);

Output.Data = reshape(Output.Data,[Settings.fov_overgrid*Output.RecoPar.DataSize(1:2) Output.Par.DataSize(3:end-1)]);

bla = ([size(Output.Data,1) size(Output.Data,2)] - Output.RecoPar.DataSize(1:2))/2+1;
Output.Data = Output.Data(bla(1):bla(1)+Output.RecoPar.DataSize(1)-1,bla(1):bla(2)+Output.RecoPar.DataSize(2)-1,:,:,:,:);

% Reconstruct noise-data non-iteratively. Cannot reconstruct noise iteratively, or could we? We actually would want to perform exactly the same reconstruction
% as for the normal data...
if(isfield(Output,'NoiseData'))
    Output.NoiseData = AOT(Output.NoiseData);
end


%% Perform Reconstruction in Slice and z-dimension

% For now just reshape them. We dont have slices or 3D-measurements for now...
Size = size(Output.Data);
Output.Data = reshape(Output.Data, [Size(1:2) prod(Size(3:4)) Size(5:end)]);
if(isfield(Output,'NoiseData'))
    Output.NoiseData = reshape(Output.NoiseData, [Size(1:2) prod(Size(3:4)) Size(5:end)]);
end


%% Conj at End

if(Settings.ConjIniSpace_flag)
    Output.Data = conj(Output.Data);
    if(isfield(Output,'NoiseData'))
        Output.NoiseData = conj(Output.NoiseData);
    end
end

%% Flip left right

% Output.Data = flip(Output.Data,2);


%% Create & Adapt Parameters


Output.RecoPar.DataSize = size(Output.Data);
Output.RecoPar.total_channel_no_reco = 1;      % I actually should make an if-condition to determine if coil combination was really done...



%% Postparations

Output = supp_UpdateRecoSteps(Output,Settings);



