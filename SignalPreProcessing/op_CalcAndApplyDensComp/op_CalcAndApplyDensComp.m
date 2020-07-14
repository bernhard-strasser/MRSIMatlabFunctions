function [MRStruct,AdditionalOut] = op_CalcAndApplyDensComp(MRStruct,NUFTOperator,Settings)
%
% op_ReadAndRecoBorjanSpiralData Read and reconstruct data from Borjan Gagoski's Spiral MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in Spiral MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         MRStruct            ...  MRStruct-struct with fields (only InTraj.GM and RecoPar have to be present)
%                                           RecoSteps:  All the settings of the performed reconstruction steps
%                                           Par:        The original parameters of the read in data
%                                           RecoPar:    The parameters after reconstruction
%                                           Data:       The non-cartesian data data
%                                           NoiseData:  If available, noise with the same scale and size as the Data is produced. Useful for SNR calculations.
%                                           OutTraj:    The Cartesian trajectory which would correspond to the Fourier transform of the Cartesian output image
%                                                       (if we did k-space gridding, those would be the k-space points we would grid to).
%                                           InTraj:     The (spiral) trajectory with which the data were measured.
% -         NUFTOperator    ...  The operator for going from an image to the NonCartesian k-Space data. Used for the AutoScale option.
% -         Settings        ...  Struct with fields
%                                           Method:         String specifying which method for calculating the DCF should be used.
%                                                           So far, only 'SpiralHoge1997AbruptChanges' implemented.
%                                           AutoScale_flag: When true, the NUFTOperator is used to construct and reconstruct an image of ones. The ground truth
%                                                           and reconstruction are compared in scale, and the ratio is applied to the DCF.
%                                                           After that, the DCF should provide properly scaled images.
%                                           Normalize_flag: If true, normalize DCF so that norm(DCF) = sqrt(numel(DCF)).
%                                           InvertDCF_flag: When true, the 1/DCF instead of DCF is computed and applied. This is rarely necessary, but
%                                                           If you want to see the k-space density, that's how you could get it.
%
% Output:
% -         MRStruct        ...  Same as for input
% -         AdditionalOut   ...  Struct with additional output, e.g. Density Compensation Function, etc.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ???

% Further remarks:



%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end

if(~isfield(Settings,'Method'))
   Settings.Method = 'SpiralHoge1997AbruptChanges';    
end
if(~isfield(Settings,'AutoScale_flag'))
   Settings.AutoScale_flag = true;    
end
if(~isfield(Settings,'Normalize_flag'))
   Settings.Normalize_flag = false;    
end
if(~isfield(Settings,'InvertDCF_flag'))
   Settings.InvertDCF_flag = false;    
end

%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(strcmpi(Settings.Method,'SpiralHoge1997AbruptChanges'))
    v1 = MRStruct.InTraj.GM;
    DCFPreG = zeros([size(v1,2) size(v1,3)]);
    for SpirPts = 2:size(v1,2)
        DCFPreG(SpirPts,:) = sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) .* ...
        abs( sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) - sqrt( v1(1,SpirPts-1,:).^2 + v1(2,SpirPts-1,:).^2 ) );
    end
    DCFPreG(isnan(DCFPreG)) = 0;
elseif(strcmpi(Settings.Method,'ConcentricRingTrajectory_Theoretical'))
    nc = MRStruct.RecoPar.nAngInts;
    R = sqrt(squeeze(MRStruct.InTraj.GM(1,1,:).^2 + MRStruct.InTraj.GM(2,1,:).^2));
    
    DCFPreG = zeros([1 nc]);
    for i=2:nc-1
        DCFPreG(i)=((R(i)/max(R)+R(i+1)/max(R))^2/4-(R(i-1)/max(R)+R(i)/max(R))^2/4)^(1)*pi;
    end
    DCFPreG(1)=(0.25*pi*(R(1)/max(R)+R(2)/max(R))^2);    
    DCFPreG(nc)=pi*((1.5*R(nc)/max(R)-0.5*R(nc-1)/max(R))^2-(0.5*R(nc)/max(R)+0.5*R(nc-1)/max(R))^2);
    DCFPreG = repmat(DCFPreG,[MRStruct.RecoPar.TrajPts 1]);
    clear nc R;
else
   st = dbstack;
   fprintf('\nWarning in %s: Did not recognize method for calculating density compensation function ''%s''.\nUse DCF = 1 instead.\n',st(1).name,Settings.Method)
   DCFPreG = ones(size_MultiDims(MRStruct.InTraj.GM,[2 3]));
end
    
%% Scale DCF
    
if(Settings.Normalize_flag)
    Scale = norm(DCFPreG(:))/sqrt(numel(DCFPreG));
elseif(Settings.AutoScale_flag)
    OnesData = ones(MRStruct.RecoPar.DataSize(1:2));
    OutOnesData = abs(NUFTOperator'*(DCFPreG(:) .* (NUFTOperator*OnesData(:)))*size(MRStruct.OutTraj.GM(:,:),2));
    OutOnesData(OutOnesData == 0) = NaN;
    Scale = nanmean(OutOnesData);
else
    %         FudgeFactor = 1.2743;     % For old trajectory
%         FudgeFactor = 0.00051078;         % For new trajectory
    FudgeFactor = 1.9634;
    Scale = max(DCFPreG(:))*2*MRStruct.RecoPar.fov_overgrid^2/FudgeFactor;
    % I dont know what these factors are. The 2*SpSpice.SimPar.fov_overgrid^2 I guessed. The FudgeFactor I got by inputting a image of ones
    % and seeing how it was scaled...
    
end
DCFPreG = DCFPreG/Scale;



%% Invert DCF if necessarry

if(Settings.InvertDCF_flag)
    DCFPreG = 1./DCFPreG; DCFPreG(isnan(DCFPreG) | isinf(DCFPreG)) = 0;
end

%% Apply DCF

if(isfield(MRStruct,'Data'))
    MRStruct.Data = MRStruct.Data .* myrepmat(single(DCFPreG(:)),size(MRStruct.Data));
end
if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData) > 1)
    MRStruct.NoiseData = MRStruct.NoiseData .* myrepmat(DCFPreG(:),size(MRStruct.NoiseData));
end
    

%% Postparations

if(nargout > 1 && exist('DCFPreG','var'))
    AdditionalOut.DCFPreG = DCFPreG;
end

% MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);


