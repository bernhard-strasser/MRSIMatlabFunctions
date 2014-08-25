function [Noise_std,Signal_ApodSmooth,Signal_SubBas] = noise_computation_new_0_3(Signal,NoiseRegion,ApodizeOrSmooth,ApodOrSmoothValue)
%
% Compute the noise of a signal (std of the signal). 
% - Do a linear fit of the sugnal and subtract that from signal. This is done because sometimes there is just 
%   simple linear baseline in the signal (eg for spectra if the water suppression did not work properly).
% - Divide the signal - LinearFit in several subregions, compute the stdev of all subregions and the stdev of the 
%   whole region.
% - Compute the ratio of the global and the local stddev. If this ratio is close to 1, there is just noise in the 
%   signal. Depending on this ratio, subtract a apodized or smoothed signal from the original signal. The 
%   smoothing/apodization is the stronger the closer the ratio is to 1, because if there is just noise inside
%   one only has to subtract something linear but no peaks.
%   With this subtraction the baseline of the signal is subtracted.
% - Compute the stddev of the Signal from which the smoothed/apodized signal was subtracted (so the baseline
%   corrected signal)
% Apodization: Filtering the fourier transformed signal (eg in the time domain) with an exponential function.

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    

% 0.3 Preparations

pause on



%% 1. Linear Fit & Subtraction


Polyfit_Signal_Real = polyfit(NoiseRegion(1):NoiseRegion(2),real(Signal(NoiseRegion(1):NoiseRegion(2))),1);
Polyfit_Signal_Imag = polyfit(NoiseRegion(1):NoiseRegion(2),imag(Signal(NoiseRegion(1):NoiseRegion(2))),1);

Signal_SubLinear = Signal(NoiseRegion(1):NoiseRegion(2)) - Polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)) - 1i*Polyval(Polyfit_Signal_Imag,NoiseRegion(1):NoiseRegion(2));



%% 1. Compute the Apodization & Smoothing Values

% Compute the mean of the Noise Signal: If it is more or less zero --> It's just noise (assuming that the range is also low, see below) --> no Apodization or Smoothing
% Noise_Mean = mean(real(Signal_SubLinear(NoiseRegion(1):NoiseRegion(2))));


% % Compute the Range (Max - Min) of the Noise. If it is low compared to the maximum of the signal --> It's just noise (assuming that Noise_mean is zero) --> no Apodization or Smoothing.
% Noise_Range = max(real(Signal(NoiseRegion(1):NoiseRegion(2)))) - min(real(Signal(NoiseRegion(1):NoiseRegion(2))));
% Noise_Range_Normalized = Noise_Range/max(real(Signal));


% Compute the Ratio of the Std of the Whole Noise Region and Subregions of the Region
TotalPeakregionPoints = NoiseRegion(2) - NoiseRegion(1) + 1;
PointsPerSubregion = 20;
TotalSubregions = floor(TotalPeakregionPoints / PointsPerSubregion);
Noise_Subregions_Std = zeros([1 TotalSubregions]);

for Subregion = 1:TotalSubregions-1
    Noise_Subregions_Std(Subregion) = std(real(Signal_SubLinear(1 + (Subregion-1)*PointsPerSubregion : Subregion*PointsPerSubregion)));
end
Noise_Subregions_Std(Subregion + 1) = std(real(Signal_SubLinear(Subregion*PointsPerSubregion : end)));
Noise_StdLocal = mean(Noise_Subregions_Std);
Noise_StdGlobal = std(real(Signal_SubLinear));
Noise_StdGlobalLocalRatio = Noise_StdGlobal / Noise_StdLocal;



if(strcmpi(ApodizeOrSmooth,'Apodize'))
    ApodOrSmoothValue = 13/sqrt(Noise_StdGlobalLocalRatio - 1)
else
    ApodOrSmoothValue = round(13/sqrt(Noise_StdGlobalLocalRatio - 1)) + 10
end



%% 2. Smooth the data or apodize

if(strcmpi(ApodizeOrSmooth,'Apodize'))
    Signal_Time = ifft(ifftshift(Signal));
    Signal_Time_Apod = Signal_Time .* exp(-ApodOrSmoothValue*(1:numel(Signal_Time))/numel(Signal_Time));
    Signal_ApodSmooth = fftshift(fft(Signal_Time_Apod));
else
    Signal_ApodSmooth_Real = transpose(smooth(real(Signal_SubLinear),ApodOrSmoothValue));
    Signal_ApodSmooth_Imag = transpose(smooth(imag(Signal_SubLinear),ApodOrSmoothValue));
    Signal_ApodSmooth = Signal;
    Signal_ApodSmooth(NoiseRegion(1):NoiseRegion(2)) = Signal_ApodSmooth_Real + 1i*Signal_ApodSmooth_Imag + ...                                                                         % Smoothing
                                                       Polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)) + 1i*Polyval(Polyfit_Signal_Imag,NoiseRegion(1):NoiseRegion(2));      % FitLinear
end




%% 3. Subtract the Smoothed or Apodized Data

Signal_SubBas = Signal - Signal_ApodSmooth;


%% 4. Compute the noise of the SubBas data

Noise_std = std(horzcat([real(Signal_SubBas(NoiseRegion(1):NoiseRegion(2))),imag(Signal_SubBas(NoiseRegion(1):NoiseRegion(2)))]));



%% DEBUG MODE: plot
% aaa = figure;
% plot(real(Signal(NoiseRegion(1):NoiseRegion(2))))
% hold on
% plot(Polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)),'r')
% hold off
% bbb = figure;
% if(strcmpi(ApodizeOrSmooth,'Apodize'))
%     plot(real(Signal(NoiseRegion(1):NoiseRegion(2))))
% else
%     plot(real(Signal_SubLinear))
% end
% hold on
% plot(real(Signal_ApodSmooth(NoiseRegion(1):NoiseRegion(2))),'g')
% hold off
% ccc = figure;
% plot(real(Signal_SubBas(NoiseRegion(1):NoiseRegion(2))))
% 
% pause
% waitforbuttonpress
% close(aaa,bbb,ccc);



%% END

pause off







