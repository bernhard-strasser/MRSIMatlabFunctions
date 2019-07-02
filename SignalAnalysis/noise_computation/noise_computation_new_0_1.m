function [Noise_std,Signal_ApodSmooth,Signal_SubBas] = noise_computation_new_0_1(Signal,NoiseRegion,ApodizeOrSmooth,ApodOrSmoothValue)
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





%% 1. Smooth the data or apodize

if(strcmpi(ApodizeOrSmooth,'Apodize'))
    Signal_Time = ifft(ifftshift(Signal));
    Signal_Time_Apod = Signal_Time .* exp(-ApodOrSmoothValue*(1:numel(Signal_Time))/numel(Signal_Time));
    Signal_ApodSmooth = fftshift(fft(Signal_Time_Apod));
    
else
    Signal_ApodSmooth_Real = transpose(smooth(real(Signal),ApodOrSmoothValue));
    Signal_ApodSmooth_Imag = transpose(smooth(imag(Signal),ApodOrSmoothValue));
    Signal_ApodSmooth = Signal_ApodSmooth_Real + 1i*Signal_ApodSmooth_Imag;
end


%% 2. Subtract the Smoothed or Apodized data

Signal_SubBas = Signal - Signal_ApodSmooth;



%% 3. Compute the noise of the SubBas data

Noise_std = std(horzcat([real(Signal_SubBas(NoiseRegion(1):NoiseRegion(2))),imag(Signal_SubBas(NoiseRegion(1):NoiseRegion(2)))]));




%% END

pause off







