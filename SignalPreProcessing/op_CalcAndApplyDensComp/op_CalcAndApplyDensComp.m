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
    if(isfield(MRStruct.RecoPar,'AssumedSequence'))
        if(strcmpi(MRStruct.RecoPar.AssumedSequence,'ViennaCRT'))
            Settings.Method = 'ConcentricRingTrajectory_Theoretical';    
        elseif(strcmpi(MRStruct.RecoPar.AssumedSequence,'BorjanSpiral'))
           Settings.Method = 'SpiralHoge1997AbruptChanges';    
        else
           Settings.Method = 'AntoinesVoronoi';                
        end
    end
end
if(~isfield(Settings,'AutoScale_flag'))
   Settings.AutoScale_flag = false;    
end
if(~isfield(Settings,'Normalize_flag'))
   Settings.Normalize_flag = false;    
end
if(~isfield(Settings,'InvertDCF_flag'))
   Settings.InvertDCF_flag = false;    
end
if(~isfield(Settings,'ApplyHammingFilter_flag'))
    Settings.ApplyHammingFilter_flag = false;
end


%%
InputName = inputname(1);
if(~isempty(InputName))
    evalin('caller',['clear ' InputName])
end

%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

DCFPreG = cell([1 numel(MRStruct.InTraj.GM)]);
if(strcmpi(Settings.Method,'SpiralHoge1997AbruptChanges'))
    for CurPart = 1:numel(MRStruct.InTraj.GM)
        v1 = MRStruct.InTraj.GM{CurPart}(:,:);
        DCFPreG{CurPart} = zeros([sum(MRStruct.RecoPar.TrajPts) 1]);
        CurPt = 0;
        for CurAngInt = 1:MRStruct.RecoPar.nAngInts
            CurPt = CurPt + 1;  % To skip the first point of the SpirPts and leave it zero
            for SpirPts = 2:MRStruct.RecoPar.TrajPts(CurAngInt)
                CurPt = CurPt + 1;
                DCFPreG(CurPt) = sqrt( v1{CurAngInt}(1,SpirPts).^2 + v1{CurAngInt}(2,SpirPts).^2 ) .* ...
                abs( sqrt( v1{CurAngInt}(1,SpirPts).^2 + v1{CurAngInt}(2,SpirPts).^2 ) - sqrt( v1{CurAngInt}(1,SpirPts-1).^2 + v1{CurAngInt}(2,SpirPts-1).^2 ) );
            end
        end
        DCFPreG{CurPart}(isnan(DCFPreG)) = 0;
    end
elseif(strcmpi(Settings.Method,'ConcentricRingTrajectory_Theoretical'))
    
    
    for CurPart = 1:numel(MRStruct.InTraj.GM)
        nc = MRStruct.RecoPar.nAngInts;
%         R = MRStruct.InTraj.GM{CurPart}(:,1,:);
%         R = sqrt(R(1,:).^2 + R(2,:).^2);
        R = MRStruct.Par.kSpaceRadii;

        DCFPreG2 = zeros([1 nc]);
        for ii=2:nc-1
            DCFPreG2(ii)=((R(ii)/max(R)+R(ii+1)/max(R))^2/4-(R(ii-1)/max(R)+R(ii)/max(R))^2/4)^(1)*pi;
        end
        DCFPreG2(1)=(0.25*pi*(R(1)/max(R)+R(2)/max(R))^2);    
        DCFPreG2(nc)=pi*((1.5*R(nc)/max(R)-0.5*R(nc-1)/max(R))^2-(0.5*R(nc)/max(R)+0.5*R(nc-1)/max(R))^2);
        DCFPreG{CurPart} = zeros([sum(MRStruct.RecoPar.TrajPts) 1]);
        CurPt = 1;
        for ii=1:nc
            DCFPreG{CurPart}(CurPt:CurPt+MRStruct.RecoPar.TrajPts(ii)-1) = repmat(DCFPreG2(ii),[MRStruct.RecoPar.TrajPts(ii) 1]);

            % Since we measure fewer points for the inner circle, need to compensate for that. Although our SNR/point increases in that case, but the scanner ADC seems
            % to adjust the scaling of the measured data depending on the ADC-dt. It seems that if you measure with longer ADC-dt, the signal does not change.
            % So effectively, with fewer points per circle, we have fewer points with the same scale as we have for more points, and thus the inner circles are
            % down-weighted.
            DCFPreG{CurPart}(CurPt:CurPt+MRStruct.RecoPar.TrajPts(ii)-1) = DCFPreG{CurPart}(CurPt:CurPt+MRStruct.RecoPar.TrajPts(ii)-1) * (max(MRStruct.RecoPar.TrajPts)./MRStruct.RecoPar.TrajPts(ii));

            CurPt = CurPt + MRStruct.RecoPar.TrajPts(ii);
        end
        clear nc R;
    end
    
elseif(strcmpi(Settings.Method,'AntoinesVoronoi'))
    for CurPart = 1:numel(MRStruct.InTraj.GM)
        DCFPreG{CurPart} = SamplingDensityVoronoi(transpose(cat(2,MRStruct.InTraj.GM{CurPart}(:,:))));
    end
    
else
   st = dbstack;
   fprintf('\nWarning in %s: Did not recognize method for calculating density compensation function ''%s''.\nUse DCF = 1 instead.\n',st(1).name,Settings.Method)
   DCFPreG = ones([sum(MRStruct.RecoPar.TrajPts) 1]);
end
    
%% Scale DCF
    

for CurPart = 1:numel(MRStruct.InTraj.GM)

    if(Settings.Normalize_flag)
        %Scale = norm(DCFPreG(:))/sqrt(numel(DCFPreG));
        Scale = norm(DCFPreG{CurPart}(:))/2300 / sqrt(numel(DCFPreG{CurPart}))*sqrt(7418); 
    elseif(Settings.AutoScale_flag)
        OnesData = ones(MRStruct.RecoPar.DataSize(1:2));
        OutOnesData = abs(NUFTOperator{CurPart}'*(DCFPreG{CurPart}(:) .* (NUFTOperator{CurPart}*OnesData(:)))*size(MRStruct.OutTraj.GM(:,:),2));
        OutOnesData(OutOnesData == 0) = NaN;
        Scale = nanmean_own(OutOnesData)*20/3;

    else
        %         FudgeFactor = 1.2743;     % For old trajectory
    %         FudgeFactor = 0.00051078;         % For new trajectory
    %     FudgeFactor = 1.9634;
    %     Scale = max(DCFPreG(:))*2*MRStruct.RecoPar.fov_overgrid^2/FudgeFactor;
    %     Scale = 64E-3;
        Scale = 1;
    %     Scale = 1;
    %     Scale = 1E-3;

        % I dont know what these factors are. The 2*SpSpice.SimPar.fov_overgrid^2 I guessed. The FudgeFactor I got by inputting a image of ones
        % and seeing how it was scaled...

    end
    DCFPreG{CurPart} = DCFPreG{CurPart}/Scale;
end

%% Apply Hamming Filter to DCF

if(Settings.ApplyHammingFilter_flag)
%     HAM=0.54+0.46*cos(2*pi*sqrt(Traj(1,:).^2+Traj(2,:).^2));
    for CurPart = 1:numel(MRStruct.InTraj.GM)
       HAM{CurPart} = transpose(0.54+0.46*cos(2*pi*sqrt(MRStruct.InTraj.GM{CurPart}(1,:).^2+MRStruct.InTraj.GM{CurPart}(2,:).^2)));
    end
    DCFPreG = cellfun(@times,DCFPreG, HAM,'uni',0);
    clear HAM;
end



%% Invert DCF if necessarry

if(Settings.InvertDCF_flag)
    for CurPart = 1:numel(MRStruct.InTraj.GM)    
        DCFPreG{CurPart} = 1./DCFPreG{CurPart}; DCFPreG{CurPart}(isnan(DCFPreG{CurPart}) | isinf(DCFPreG{CurPart})) = 0;
    end
end

%% Apply DCF

for CurPartEnc = 1:MRStruct.RecoPar.nPartEnc
    
    if(numel(MRStruct.Data) ~= numel(DCFPreG))      % If for each partition we have identical trajectories, just maybe fewer
        DCFPart = 1;
    else                                            % If each partition has unique trajectories
        DCFPart = CurPartEnc;
    end
    
    if(isfield(MRStruct,'Data'))
            MRStruct.Data{CurPartEnc} = MRStruct.Data{CurPartEnc} .* single(DCFPreG{DCFPart}(1:size(MRStruct.Data{CurPartEnc},1)));
    end
    if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData{1}) > 1)
            MRStruct.NoiseData{CurPartEnc} = MRStruct.NoiseData{CurPartEnc} .* single(DCFPreG{DCFPart}(1:size(MRStruct.NoiseData{CurPartEnc},1)));
    end

    
    
    
end  

%% Postparations

if(nargout > 1 && exist('DCFPreG','var'))
    AdditionalOut.DCFPreG = DCFPreG;
end

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);


