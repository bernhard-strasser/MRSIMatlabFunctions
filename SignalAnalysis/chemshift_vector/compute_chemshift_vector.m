function [chemshift_vector, FIDTime, freq_vector, bandwidth_frequency, step_frequency] = compute_chemshift_vector(LarmorFreqOrPar,dwelltime,vecSize,CenterAroundPPM)
% In units of:
% water_frequency: Hz
% dwelltime: s
% vecSize: -


%% Preparations

if(~exist('CenterAroundPPM','var'))
    CenterAroundPPM = 4.65;
end
if(isstruct(LarmorFreqOrPar))
    
    if(isfield(LarmorFreqOrPar,'RecoPar'))
        LarmorFreqOrPar = LarmorFreqOrPar.RecoPar;
    end
    if(isfield(LarmorFreqOrPar,'Par'))
        LarmorFreqOrPar = LarmorFreqOrPar.Par;
    end    
    LarmorFreq = LarmorFreqOrPar.LarmorFreq;
    if(~exist('dwelltime','var') || isempty(dwelltime))
        dwelltime = LarmorFreqOrPar.Dwelltimes(1)/1E9;
    end
    if(~exist('vecSize','var') || isempty(vecSize))
        vecSize = LarmorFreqOrPar.vecSize(1);
    end
else
    LarmorFreq = LarmorFreqOrPar;
end

%% Compute chemshift_vector

% 1) Calculate the spectral bandwidth = 1/dwelltime.
% 2) The stepsize in frequency is then bandwidth/vecsize.
% 3) And so the frequency vector goes from (vecSize/2-1)*step_frequency to -vecSize/2*step_frequency in steps of -step_frequency.
%    However, need ceil and floor in case vecSize is odd. 
%    Above formulae are actually a definition for even vecSizes. Example: vecSize = 512, aboce formula gives (255:-1:-256)*step_frequency. We could also define it
%    as (256:-1:-255)*step_frequency, or (255.5:-1:-255.5)*step_frequency. The last definition would be centered around 0, but the first definition is similar to the
%    k-space definition, where for even matrix sizes, k-space is also not centered around 0,  but is (-N/2):1:(N/2-1).
% 4) Need to convert the frequency into chemical shift, and add CenterAroundPPM, if we acquire our spectra e.g. centered around the water-peak, i.e. 4.65 ppm.


bandwidth_frequency = 1/dwelltime;
step_frequency = bandwidth_frequency / vecSize;
freq_vector = ((ceil(vecSize/2)-1):-1:-floor(vecSize/2)) * step_frequency;

chemshift_vector = 10^6*freq_vector / LarmorFreq + CenterAroundPPM;


%% FID Time

FIDTime = (0:dwelltime:(dwelltime*(vecSize-1)));

