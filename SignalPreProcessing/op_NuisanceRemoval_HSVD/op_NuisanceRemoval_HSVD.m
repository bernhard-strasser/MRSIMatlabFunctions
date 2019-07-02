function [csi, csi_Nuis] = op_NuisanceRemoval_HSVD(csi, Settings,mask,quiet_flag)
%
% Remove nuisance (fat, water) signal from spectroscopy data using HSVD
%
% This function was written by Bernhard Strasser, October 2018.
% The functions used within this function are written by Chao Ma.
%
% The function ...
%
%
% [csi, csi_Nuis] = NuisanceRemoval_HSVD(csi, Settings,mask,quiet_flag)
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
% File dependancy:

% Remark: If the data were conjugated, the chemical shift definition is the negative of the conventional one:
% Delta_conj := -Delta := (f_TMS - f)/f_TMS * 10^6. So positive frequencies starting from the water-frequency (which is set to 0
% because we look at f - f_Water) are towards lower chemical shifts, and negative frequencies towards higher chemical shifts. That
% means that NAA e.g. has a POSITIVE frequency (normally it would have negative).
% Set the flag 'Settings.ConjugationFlag' apropriately!


%% Preparations

if(~isfield(Settings,'Debug_flag'))
   Settings.Debug_flag = false; 
end

% ConjugationSign: If data were conjugated (to flip spectrum), the calculation of the frequency from the ppm is different:
% it is f_conj = LarmorFreq/10^6 * (4.65 - Delta) instead of f = LarmorFreq/10^6 * (Delta - 4.65).
% By default, assume the data were conjugated
if(~isfield(Settings,'ConjugationFlag'))
    if(~isfield_recursive(csi,'RecoPar.ConjugationFlag'))
        Settings.ConjugationSign = -1;
    else
        Settings.ConjugationSign = (-1)^csi.RecoPar.ConjugationFlag;        
    end
else
    Settings.ConjugationSign = (-1)^Settings.ConjugationFlag;            
end
% Turn off the annoying rank deficiency warnings from the mldivide ( "/" operator). Turn it on again after the function ran.
warning('off','MATLAB:rankDeficientMatrix')

%% Some definitions


dim    = size(csi.Data);
dt     = csi.Par.Dwelltimes(1)*1e-9;


%% water removal 

Factor = csi.Par.LarmorFreq * 1e-6;

% Signal selective HSVD params
selParams               = struct('dt',dt,'n',Settings.NoOfSingVals); % n: Number of singular values that should be kept in the truncated SVD.
optWat                  = struct('fSel2',Settings.ConjugationSign*(Settings.WaterPPMs - 4.65)*Factor,'maxT2',Settings.WaterT2s); % Chaos chemical shift definition
optLip                  = struct('fSel2',Settings.ConjugationSign*(Settings.LipidPPMs - 4.65)*Factor,'maxT2',Settings.LipidT2s); % is negative of normal!!!
optMeta.fSel2           = Settings.ConjugationSign*(Settings.MetaboPPMs - 4.65)*Factor;                                     % Parameters of the metabolites.
if(~isempty(Settings.OtherPPMs))
    optOther = struct('fSel2',Settings.ConjugationSign*(Settings.OtherPPMs - 4.65)*Factor,'maxT2',Settings.OtherT2s);
else
    optOther                = [];    % If voxel is not in mask, process data with that settings. If empty, skip voxel.
end
selParams.signalType    = 'water';
selParams.NtEcho        = 0;
if(isfield(Settings,'realT2'))
    selParams.realT2 = Settings.realT2;
end

if(~quiet_flag)
    fprintf('\n\nRemove nuisance signals . . . ');
end

if(sum(sum(sum(mask))) > 200 && license('test','Parallel_Computing_Toolbox'))
    Parallel_flag = true;
else
    Parallel_flag = [];
end
    
SizeCsi = size(csi.Data);
csi.Data = csi.Data(:,:,:,:,:);       % Put all dimensions > 5 to the 5th dimension
NuisTic = tic;
if(nargout > 1)
    csi_Nuis = zeros(size(csi.Data));
end
if(Settings.Debug_flag)
    bak = csi.Data;
end
for OtherDimInd = 1:prod(SizeCsi(5:end))
    for slc = 1:size(csi.Data,3)
        if(nargout > 1 || Settings.Debug_flag)
            [csi.Data(:,:,slc,:,OtherDimInd),csi_Nuis(:,:,slc,:,OtherDimInd)] = nsRm_bstrchanged ...
            (squeeze(csi.Data(:,:,slc,:,OtherDimInd)), mask(:,:,slc), zeros(dim(1:2)), selParams,optWat, optLip, optMeta, optOther,Parallel_flag);
        else
            csi.Data(:,:,slc,:,OtherDimInd) = nsRm_bstrchanged ...
            (squeeze(csi.Data(:,:,slc,:,OtherDimInd)), mask(:,:,slc), zeros(dim(1:2)), selParams,optWat, optLip, optMeta, optOther,Parallel_flag);            
        end
    end
end
csi.Data = reshape(csi.Data,SizeCsi);
if(~quiet_flag)
    fprintf('%s seconds.',num2str(toc(NuisTic)))
end

if(Settings.Debug_flag)
    for CurVox = 1:numel(Settings.DebugVoxels)
        chemy = compute_chemshift_vector_1_2(csi.Par.LarmorFreq,csi.Par.Dwelltimes(1)/10^9,csi.Par.vecSize);
        figure; 
        plot(chemy,squeeze(abs(fftshift(fft(bak(Settings.DebugVoxels{CurVox}(1),Settings.DebugVoxels{CurVox}(2),Settings.DebugVoxels{CurVox}(3),:,1,1))))))
        hold on
        plot(chemy,squeeze(abs(fftshift(fft(csi_Nuis(Settings.DebugVoxels{CurVox}(1),Settings.DebugVoxels{CurVox}(2),Settings.DebugVoxels{CurVox}(3),:,1,1))))),'r')
        hold off
    end
end



%% Postparations

warning('on','MATLAB:rankDeficientMatrix')

