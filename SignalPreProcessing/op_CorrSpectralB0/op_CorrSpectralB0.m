function [InData, AdditionalOut] = op_CorrSpectralB0(InData,B0,Mask,Settings)
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

if(isempty(Mask))
    clear Mask
end
if(~exist('Settings','var'))
    Settings = [];
end
if(~isfield(Settings,'Debug_flag'))
    Settings.Debug_flag = false;
end
if(~isfield(Settings,'Overdiscrete_flag'))
    Settings.Overdiscrete_flag = false;
end
if(Settings.Overdiscrete_flag && ~isfield(Settings,'OverdiscreteSize'))
    Settings.OverdiscreteSize = size_MultiDims(B0.B0Map,1:3);
end
if(~Settings.Overdiscrete_flag || ~isfield(Settings,'Downsample_flag'))
    Settings.Downsample_flag = false;
end
if(~isfield(Settings,'ReverseB0_flag'))
    Settings.ReverseB0_flag = false;
end
if(~isfield(Settings,'RoundB0ToIntVecPts'))
    Settings.RoundB0ToIntVecPts = false;
end
if(~isfield(InData.Par,'DataSize'))
    InData.Par.DataSize = size(InData.Data);
end
if(~isfield(InData,'RecoPar'))
    InData.RecoPar = InData.Par;
end
if(~isfield(InData.RecoPar,'DataSize'))
    InData.RecoPar.DataSize = size(InData.Data);
end


%% Calculate B0-Map 
% From Residual water and lipids in case B0-map is not provided

if(~exist('B0','var') || isempty(B0))
    TestIn.csi = InData.Data;
    TestIn.mask = InData.Mask;
    Sett = struct('PeakSearchPPM',[4.7],'PolyfitRegion',[4.4 5.0],'PeakSearchRangePPM',0.5); 
    Sett.LarmorFreq = InData.RecoPar.LarmorFreq; Sett.Dwelltime = InData.RecoPar.Dwelltimes(1); Sett.vecsize = InData.RecoPar.vecSize; 
    TestIn.csi = InData.Data; TestIn.mask = InData.Mask;
    [TestOut,ShiftMap] = FrequencyAlignment(TestIn,Sett,4,2); 
    InData.Data = TestOut;
end


%% Reverse B0

if(Settings.ReverseB0_flag)
    B0.B0Map = -1*B0.B0Map;
end


%% DEBUG: Plot Spectra Before B0Corr

if(Settings.Debug_flag)
    bla = InData.Data(:,:,1,:,1);
    bla = bla .* myrepmat(InData.Mask,size(bla));
    bla = abs(fftshift(fft(bla,[],4),4));
    bla = permute(bla,[4 1 2 3]);
    chemy = compute_chemshift_vector_1_2(InData.RecoPar.LarmorFreq,InData.RecoPar.Dwelltimes(1)/10^9,InData.RecoPar.vecSize);
    figure; plot(chemy,bla(:,:))
    figure; plot(chemy,sum(bla(:,:),2))
    clear bla chemy
end

%% Get to Same Resolution for B0-map and InData

if(Settings.Overdiscrete_flag)
    InData.Data = ZerofillOrCutkSpace(InData.Data,[Settings.OverdiscreteSize InData.RecoPar.DataSize(4:end)],1);
    if(isfield(InData,'NoiseData'))
        InData.NoiseData = ZerofillOrCutkSpace(InData.NoiseData,[Settings.OverdiscreteSize InData.RecoPar.DataSize(4:end)],1);       
    end
end

if(ndims(B0.B0Map) == 3)
    CurB0 = imresize3(B0.B0Map,size_MultiDims(InData.Data,1:3));
else   
    CurB0 = imresize(B0.B0Map,size_MultiDims(InData.Data,1:2));
end
if(exist('Mask','var'))
    if(ndims(Mask) == 3)
        Mask = imresize3(double(Mask),size(CurB0),'nearest');
    else
        Mask = imresize(double(Mask),size(CurB0),'nearest');        
    end
    CurB0 = CurB0 .* Mask;
end

% CAREFUL: FROM HERE ON InData.RecoPar.DataSize MIGHT HAVE NOT THE REAL SIZE OF THE InData.Data !!!


%% Apply B0Mat to D1Cart
% Round to shift only integer number of points
HzPerPt = 10^9/InData.RecoPar.Dwelltimes(1) / InData.RecoPar.vecSize;
if(Settings.RoundB0ToIntVecPts)
    CurB0 = round(CurB0/HzPerPt)*HzPerPt;
end

% Remove NaN's
CurB0(isnan(CurB0)) = 0;

AdditionalOut.B0Map = CurB0;
time   = (0:InData.RecoPar.DataSize(4)-1)*InData.RecoPar.Dwelltimes(1)/10^9;
AdditionalOut.TimeVec = time;
B0CorrMat_Spec = exp(myrepmat(2*pi*1i*CurB0,size(InData.Data)) .* myrepmat(time,size(InData.Data)));

InData.Data = InData.Data .* B0CorrMat_Spec;
if(isfield(InData,'NoiseData'))
    InData.NoiseData = InData.NoiseData .* B0CorrMat_Spec;    
end

%% Save B0CorrMat_Spec

if(nargout > 1)
    AdditionalOut.B0CorrMat_Spec = B0CorrMat_Spec;
end


%% Downsample Data
if(Settings.Overdiscrete_flag)
    if(Settings.Downsample_flag)
        InData.Data = ZerofillOrCutkSpace(InData.Data,InData.RecoPar.DataSize,1);
        if(isfield(InData,'NoiseData'))
            InData.NoiseData = ZerofillOrCutkSpace(InData.NoiseData,InData.RecoPar.DataSize,1);       
        end
    else
        InData.RecoPar.DataSize(1:3) = size_MultiDims(B0.B0Map,1:3);
    end
end

%% DEBUG: Plot Spectra After B0Corr

if(Settings.Debug_flag)
    plot_SpecOfAllVoxels(InData);
end


%% Postparations

InData = supp_UpdateRecoSteps(InData,Settings);

