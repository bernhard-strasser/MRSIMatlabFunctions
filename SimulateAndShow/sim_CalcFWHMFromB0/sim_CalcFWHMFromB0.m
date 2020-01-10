function [FWHMMap, AdditionalOut] = sim_CalcFWHMFromB0(B0,Settings)
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
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadSimMRSIDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadSimMRSIDataSets          ...     
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
    Settings = [];
end
if(~isfield(Settings,'Debug_flag'))
    Settings.Debug_flag = false;
end
if(~isfield(Settings,'DownsampleSize'))
    Settings.DownsampleSize = ceil([size(B0.B0Map,1) size(B0.B0Map,2) size(B0.B0Map,3)]/4);
end



%% Simulate Data

SimMRSIData.Par.Dwelltimes = 5E6;    % 100 Hz
SimMRSIData.Par.vecSize = 200;
NatLW = 2*pi;
t = 0:SimMRSIData.Par.Dwelltimes:(SimMRSIData.Par.vecSize-1)*SimMRSIData.Par.Dwelltimes;
t = t/1E9;
HzPerPt = 1E9/SimMRSIData.Par.Dwelltimes/SimMRSIData.Par.vecSize;


SimMRSIData.Data = exp(-t*NatLW);
SimMRSIData.Data = myrepmat(SimMRSIData.Data,[size(B0.B0Map,1) size(B0.B0Map,2) size(B0.B0Map,3) SimMRSIData.Par.vecSize]);


%% Apply B0-Shifts

SimMRSIData = op_CorrSpectralB0(SimMRSIData,B0,B0.BrainMask);



%% Downsample & Filter Simulated MRSI Data

SimMRSIData.Data = ZerofillOrCutkSpace(SimMRSIData.Data,[Settings.DownsampleSize SimMRSIData.Par.vecSize],1);
SimMRSIData.Data = HammingFilter(SimMRSIData.Data,[1 2 3],100,[],0);
SimMRSIData.RecoPar.DataSize = size(SimMRSIData.Data);
SimMRSIData.BrainMask = imresize3(double(B0.BrainMask),SimMRSIData.RecoPar.DataSize(1:3),'nearest');

%% Calculate FWHM

SpecData = fftshift(fft(SimMRSIData.Data,[],4),4);
[MaxMap, PeakPtMap] = max(real(SpecData),[],4);
FWHMMap = zeros(size(MaxMap));

for xx = 1:size(FWHMMap,1)
    for yy = 1:size(FWHMMap,2)
        for zz = 1:size(FWHMMap,3)
            if(~SimMRSIData.BrainMask(xx,yy,zz))
                continue;
            end
            Tmp1 = FindClosestIndex(squeeze(real(SpecData(xx,yy,zz,1:PeakPtMap(xx,yy,zz)))),MaxMap(xx,yy,zz)/2);
            Tmp2 = FindClosestIndex(squeeze(real(SpecData(xx,yy,zz,PeakPtMap(xx,yy,zz)+1:end))),MaxMap(xx,yy,zz)/2);
            HalfMaxPts = [Tmp1{1}(1) Tmp2{1}(1)+PeakPtMap(xx,yy,zz)];
            FWHMMap(xx,yy,zz) = (max(HalfMaxPts) - min(HalfMaxPts))*HzPerPt;
        end
    end
    
end



%% Postparations

AdditionalOut.SimMRSIData = SimMRSIData;



