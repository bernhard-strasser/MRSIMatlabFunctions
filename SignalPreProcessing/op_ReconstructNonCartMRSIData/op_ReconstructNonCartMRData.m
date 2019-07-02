function [Data_i, AdditionalOut] = op_ReconstructNonCartMRData(Data_k,InTrajectory,OutTrajectory,RecoPar,Settings)
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
   Settings.Phaseroll_flag = true;
   Settings.DensComp_flag = true;   
    
end
if(~isfield(Settings,'Phaseroll_flag'))
   Settings.Phaseroll_flag = true;    
end
if(~isfield(Settings,'DensComp_flag'))
   Settings.DensComp_flag = true;    
end
if(~isfield(Settings,'ConjInBegin_flag'))
   Settings.ConjInBegin_flag = false;    
end
if(~isfield(Settings,'ConjAtEnd_flag'))
   Settings.ConjAtEnd_flag = true;    
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



%% Conj in Beginning

if(Settings.ConjInBegin_flag)
    Data_k = conj(Data_k);
end

%% Try to correct for "tilted trajectory" by doing phaseroll

% Save the data for reconstructing the pseudo-pcg case
if(Settings.Phaseroll_flag)

    nTI = RecoPar.nTempInt;
    vs = RecoPar.vecSize;
    ns = RecoPar.TrajPts;
    nc = RecoPar.nAngInts;
    nrew = RecoPar.RewPts;
    ncha = size(Data_k,6);

    timeoffset = 0:(ns-1);
    timeoffset = repmat(transpose(timeoffset),[1 vs]);
    Freq = ((0:vs-1)/vs-0.5)/(nrew + ns)*nTI;
    Freq = repmat(Freq,[ns 1]);    

    % This comes from:
    % timeoffset = (0:(ns-1))*RecoPar.ADC_Dt/10^6;
    % sBW = nTI/((nrew + ns)*RecoPar.ADC_Dt/10^6);
    % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
    % RecoPar.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
    % the rest is basically the same (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to ((0:vs-1)/vs-0.5), and the other constants are
    % the same anyway
    
    phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    
    phasecorr = myrepmat(phasecorr,size(Data_k));
    
    bla = Data_k(:,:,:,:,1,:,:);
    Data_k = fft(fftshift(conj(phasecorr).*fftshift(ifft(Data_k,[],5),5),5),[],5);
    Data_k(:,:,:,:,1,:,:) = bla;

    TiltTrajMat = reshape(phasecorr(:,:,1,1,:,1),[RecoPar.TrajPts*RecoPar.nAngInts RecoPar.vecSize]);   
    
else
    
	TiltTrajMat = ones([RecoPar.TrajPts*RecoPar.nAngInts RecoPar.vecSize]);

end

if(nargout > 1)
    AdditionalOut.TiltTrajMat = TiltTrajMat;
end




%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(Settings.DensComp_flag)
    % FudgeFactor = 1.2743;     % For old trajectory
    FudgeFactor = 0.00051078;         % For new trajectory
    FudgeFactor = 1.9634;
    
    v1 = InTrajectory;
    DCFPreG = zeros([size(v1,2) size(v1,3)]);
    for SpirPts = 2:size(v1,2)
        DCFPreG(SpirPts,:) = sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) .* ...
        abs( sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) - sqrt( v1(1,SpirPts-1,:).^2 + v1(2,SpirPts-1,:).^2 ) );
    end
    DCFPreG(isnan(DCFPreG)) = 0;
    DCFPreG = DCFPreG / max(DCFPreG(:))*2*Settings.fov_overgrid^2/FudgeFactor;  %
    % I dont know what these factors are. The 2*SpSpice.SimPar.fov_overgrid^2 I guessed. The FudgeFactor I got by inputting a image of ones
    % and seeing how it was scaled...

    clear v1
    
    Data_k = Data_k .* myrepmat(DCFPreG,size(Data_k));

    
end
if(nargout > 1 && exist('DCFPreG','var'))
    AdditionalOut.DCFPreG = DCFPreG;
end


%% Calculate sft2-Operator

% sft Operator
% Reshape trajectories to expected shape

% Collapse data to a matrix (from [nAngInt x nTrajPoints x nTempInt*vecSize x nCha x nPart*nSlc] to [nAngInt*nTrajPoints x Rest])
SizeData_k = size(Data_k); SizeData_k = cat(2,SizeData_k,ones([1 5-numel(SizeData_k)]));
Data_k = reshape(Data_k,[prod(SizeData_k(1:2)) prod(SizeData_k(3:end))]);   



sft2_Oper = sft2_Operator(transpose(squeeze(OutTrajectory(:,:))*RecoPar.DataSize(1)),transpose(InTrajectory(:,:)),1);



%% Calculate B0-Correction of Spiral Data in Spatial Domain
% if(Ctrl.SimPar.B0_flag && Ctrl.SimPar.TiltedTraj_flag)
%     t   = (0:SpSpice.D2.RecoPar.DataSize(1)-1)*SpSpice.D2.RecoPar.ADC_Dt/10^6;
%     t = repmat(t,[1 1 SpSpice.D2.RecoPar.nAngInts]); t = t(:);
%     SpSpice.Reco.B0CorrMat_Spatial = exp(myrepmat(-2*pi*1i*SpSpice.D2.B0Map(:),size(sft2_Operator_ForSpir_TrajIdeal)) .* myrepmat(t,size(sft2_Operator_ForSpir_TrajIdeal)));
% 
%     % SpSpice.Spi2Cart.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;
%     % SpSpice.Spi2Cart_NoTiltCorr.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;
% 
% else
%     SpSpice.Reco.B0CorrMat_Spatial = ones(size(sft2_Operator_ForSpir_TrajIdeal));
% end
% if(Ctrl.RecoPar.Correct4B0_flag)
%     sft2_Operator_ForSpir_TrajIdeal_B0Corr = sft2_Operator_ForSpir_TrajIdeal .* SpSpice.Reco.B0CorrMat_Spatial;
% else
%     sft2_Operator_ForSpir_TrajIdeal_B0Corr = sft2_Operator_ForSpir_TrajIdeal;
% end




%% Apply sft2-Operator (Fourier-Transform from spiral k-space --> Cartesian image 

Data_i = Data_k; clear Data_k
Data_i = sft2_Oper' * Data_i * size(OutTrajectory(:,:),2);
Data_i = reshape(Data_i,[Settings.fov_overgrid*RecoPar.DataSize(1:2) SizeData_k(3:end)]);

bla = ([size(Data_i,1) size(Data_i,2)] - RecoPar.DataSize(1:2))/2+1;
Data_i = Data_i(bla(1):bla(1)+RecoPar.DataSize(1)-1,bla(1):bla(2)+RecoPar.DataSize(2)-1,:,:,:,:);

if(nargout > 1)
    AdditionalOut.sft2_Oper = sft2_Oper;
end


%% Perform Reconstruction in Slice and z-dimension

% For now just reshape them. We dont have slices or 3D-measurements for now...
Size = size(Data_i);
Data_i = reshape(Data_i, [Size(1:2) prod(Size(3:4)) Size(5:end)]); 



%% Conj at End

if(Settings.ConjAtEnd_flag)
    Data_i = conj(Data_i);
end


%% Flip left right

% Data_i = flip(Data_i,2);


