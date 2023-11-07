function [Input, AdditionalOut] = op_CutOrExtendSpectra(Input,Settings,DataForExtending)
%
% op_CutOrExtendSpectra Cuts or Extends Spectra to Certain ppm Range (Extension by Zerofilling, or Using Additional Data) 
%
% This function was written by Bernhard Strasser, March 2019.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         Input                    ...     The Input. Must have fields 'Data', and 'Par'.
% -         CoilWeightMap             ...     The Coil-weighting-map. Must have field 'Data'
% -         Settings          ...            Structure to specify how the coil combination should be performed.
%
% Output:
% -         Input                      ...     The Output (Input is overwritten). Will have field 'RecoPar' additionally.
% -         AdditionalOut                        ...  Variable for storing additional output (if asked for by user).   
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: For now, the coil index in CoilWeightMap is the first, but in Input it's the last. Change that later!
% Also, the function cannot handle 3D / multislice data for now. And the channel index is fixed for now.


% This function expects the input to be of form
% ??? 


%% 0. Preparations

if(~exist('Settings','var'))
    Settings = [];
end
if(~isfield(Settings,'PerformFFT_flag'))
    Settings.PerformFFT_flag = true;
end
if(~isfield_recursive(Input,'Par.FIDDimension'))
    Input.Par.FIDDimension = find(Input.Par.DataSize == Input.Par.vecSize);
    if(numel(Input.Par.FIDDimension) > 1)
        error('Error in op_CutOrExtendSpectra: Cannot auto-find FID dimension. Please provide manually!')
    end
end
if(~isfield(Input,'RecoPar'))
    Input.RecoPar = Input.Par;
end
if(~isfield(Input.RecoPar,'LowestPPMPoint'))
   Input.RecoPar.LowestPPMPoint = [];
end


Settings.CutOrExtendToPpm = sort(Settings.CutOrExtendToPpm); Settings.CutOrExtendToPpm = Settings.CutOrExtendToPpm(end:-1:1);
Dims = numel(Input.RecoPar.DataSize);


%% FFT from Time --> Frequency
if(Settings.PerformFFT_flag)
    Input.Data = fftshift(fft(Input.Data,[],Input.Par.FIDDimension),Input.Par.FIDDimension);
    if(isfield(Input,'NoiseData'))
        Input.NoiseData = fftshift(fft(Input.Data,[],Input.Par.FIDDimension),Input.Par.FIDDimension);
    end
end

%% Decide whether to cut or extend for each side of spectrum

ppm = compute_chemshift_vector(Input.RecoPar.LarmorFreq,Input.RecoPar.Dwelltimes(1)/10^9,Input.RecoPar.vecSize,Input.RecoPar.LowestPPMPoint);

CutFlag = [max(Settings.CutOrExtendToPpm)<max(ppm) min(Settings.CutOrExtendToPpm)>min(ppm)];

%% Left Spectral Side

if(CutFlag(1))
    % Find index on left side of spectra where we want to cut
    Ind = FindClosestIndex(ppm,Settings.CutOrExtendToPpm(1),1E-6);
    % Cut data Out = In(:,:,...,Ind:end,:,:,...); Do that quite complicated, because we don't know the size of the input array a priori!
    S.subs = repmat({':'},[1 Dims]); S.subs{Input.Par.FIDDimension} = Ind{1}:Input.RecoPar.DataSize(Input.Par.FIDDimension); S.type = '()';
    Input.Data = subsref(Input.Data,S);
    if(isfield(Input,'NoiseData'))
        Input.NoiseData = subsref(Input.NoiseData,S);
    end
else
    if(exist('DataForExtending','var'))
        ppmExt = compute_chemshift_vector(...
        DataForExtending.RecoPar.LarmorFreq,DataForExtending.RecoPar.Dwelltimes(1)/10^9,DataForExtending.RecoPar.vecSize);
        IndExt = FindClosestIndex(ppmExt,[Settings.CutOrExtendToPpm(1) max(ppm)],1E-6);
        S.subs = repmat({':'},[1 Dims]); S.subs{Input.Par.FIDDimension} = IndExt{1}(1):IndExt{2}(1)-1; S.type = '()';
        AppendData = subsref(DataForExtending.Data,S);
    else
        ExtendByPoints = round((Settings.CutOrExtendToPpm(1)-max(ppm))/abs(ppm(1)-ppm(2)));
        AppendSize = Input.RecoPar.DataSize; AppendSize(Input.Par.FIDDimension) = ExtendByPoints;
        AppendData = zeros(AppendSize);
    end
    Input.Data = cat(Input.Par.FIDDimension,AppendData,Input.Data);
end
ExtendedSize = size(Input.Data,Input.Par.FIDDimension) - Input.RecoPar.DataSize(Input.Par.FIDDimension);

%% Right Spectral Side

if(CutFlag(2))
    % Find index on left side of spectra where we want to cut
    Ind = FindClosestIndex(ppm,Settings.CutOrExtendToPpm(2),1E-6);
    % Cut data Out = In(:,:,...,Ind:end,:,:,...); Do that quite complicated, because we don't know the size of the input array a priori!
    S.subs = repmat({':'},[1 Dims]); S.subs{Input.Par.FIDDimension} = 1:(Ind{1}+ExtendedSize); S.type = '()';
    Input.Data = subsref(Input.Data,S);
    if(isfield(Input,'NoiseData'))
        Input.NoiseData = subsref(Input.NoiseData,S);
    end
    Input.RecoPar.LowestPPMPoint = ppm(Ind{1});       % Because we might cut the spectra not centered around the reference frequency
else
    if(exist('DataForExtending','var'))
        ppmExt = compute_chemshift_vector(...
        DataForExtending.RecoPar.LarmorFreq,DataForExtending.RecoPar.Dwelltimes(1)/10^9,DataForExtending.RecoPar.vecSize);
        IndExt = FindClosestIndex(ppmExt,[Settings.CutOrExtendToPpm(2) min(ppm)],1E-6);
        S.subs = repmat({':'},[1 Dims]); S.subs{Input.Par.FIDDimension} = IndExt{2}(1)+1:IndExt{1}(1); S.type = '()';
        AppendData = subsref(DataForExtending.Data,S);
    else
        ExtendByPoints = round((min(ppm)-Settings.CutOrExtendToPpm(2))/abs(ppm(1)-ppm(2)));
        AppendSize = Input.RecoPar.DataSize; AppendSize(Input.Par.FIDDimension) = ExtendByPoints;
        AppendData = zeros(AppendSize);
    end
    Input.Data = cat(Input.Par.FIDDimension,Input.Data,AppendData);
end


%% iFFT from Frequency --> Time
if(Settings.PerformFFT_flag)
    Input.Data = ifft(ifftshift(Input.Data,Input.Par.FIDDimension),[],Input.Par.FIDDimension);
end


%% Update Parameters

OrigVS = Input.RecoPar.vecSize;
Input.RecoPar.vecSize = size(Input.Data,Input.Par.FIDDimension);
Input.RecoPar.DataSize(Input.Par.FIDDimension) = Input.RecoPar.vecSize;
Input.RecoPar.Dwelltimes = Input.RecoPar.Dwelltimes * (OrigVS/Input.RecoPar.vecSize);


