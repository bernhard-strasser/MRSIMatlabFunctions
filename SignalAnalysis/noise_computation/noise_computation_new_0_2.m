function [Noise_std,Signal_ApodSmooth,Signal_SubBas] = noise_computation_new_0_2(Signal,NoiseRegion,ApodizeOrSmooth,ApodOrSmoothValue)
%
% Compute the noise of a signal (std of the signal). It divides the signal in several subregions. If there is a higher std in a certain subregion than in the others it will be excluded.
% If the standard deviations are varying very severe (if the standard deviation of all the standard deviations is high) --> abort program with error
%
%

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    

% 0.3 Preparations

pause on




%% 1. Compute the Apodization & Smoothing Values

% Compute the Ratio of the Std of the Whole Noise Region and Subregions of the Region
TotalPeakregionPoints = NoiseRegion(2) - NoiseRegion(1) + 1;
PointsPerSubregion = 20;
TotalSubregions = floor(TotalPeakregionPoints / PointsPerSubregion);
Noise_Subregions_Std = zeros([1 TotalSubregions]);

for Subregion = 1:TotalSubregions-1
    Noise_Subregions_Std(Subregion) = std(real(Signal(NoiseRegion(1) + (Subregion-1)*PointsPerSubregion : NoiseRegion(1) + Subregion*PointsPerSubregion)));
end
Noise_Subregions_Std(Subregion + 1) = std(real(Signal(NoiseRegion(1) + Subregion*PointsPerSubregion : NoiseRegion(2))));
Noise_StdLocal = mean(Noise_Subregions_Std);
Noise_StdGlobal = std(real(Signal(NoiseRegion(1):NoiseRegion(2))));
Noise_StdGlobalLocalRatio = Noise_StdGlobal / Noise_StdLocal



if(strcmpi(ApodizeOrSmooth,'Apodize'))
    ApodOrSmoothValue = 15/sqrt(Noise_StdGlobalLocalRatio - 1)
else
    ApodOrSmoothValue = round(20/sqrt(Noise_StdGlobalLocalRatio - 1)) + 14
end



%% 2. Smooth the data or apodize

if(strcmpi(ApodizeOrSmooth,'Apodize'))
    Signal_Time = ifft(ifftshift(Signal));
    Signal_Time_Apod = Signal_Time .* exp(-ApodOrSmoothValue*(1:numel(Signal_Time))/numel(Signal_Time));
    Signal_ApodSmooth = fftshift(fft(Signal_Time_Apod));
    
else
    Signal_ApodSmooth_Real = transpose(smooth(real(Signal),ApodOrSmoothValue));
    Signal_ApodSmooth_Imag = transpose(smooth(imag(Signal),ApodOrSmoothValue));
    Signal_ApodSmooth = Signal_ApodSmooth_Real + 1i*Signal_ApodSmooth_Imag;
end


%% 3. Subtract the Smoothed or Apodized data

Signal_SubBas = Signal - Signal_ApodSmooth;



%% 4. Compute the noise of the SubBas data

Noise_std = std(horzcat([real(Signal_SubBas(NoiseRegion(1):NoiseRegion(2))),imag(Signal_SubBas(NoiseRegion(1):NoiseRegion(2)))]));




%% END

pause off







