function [FID,Spectrum] = Simulate_FID_Spectra_1_3(Chemshift,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,LarmorFreq,SmoothFIDSpan)
% Simulates FIDs. Input: Angular Frequency [rad/s], Zero-Order-Phase [Degree], AcquisitionDelay[s], T2[s], Magnitude of Signal at t = 0, SNR, dwelltime[s], vecSize, Span to Smooth FID

%% Comments, Definitions, Initializations
% Define Variables

%water_frequency = 297*10^6; % @ 7T
%water_frequency = LarmorFreq * (1 + 4.65/10^6) * 2*pi;
omega = LarmorFreq * (1 + (Chemshift - 4.65)/10^6) * 2*pi;

% Define time
t_end = AcqDelay + dwelltime*(vecSize - 1);
t=AcqDelay:dwelltime:t_end;

% Define Standard Deviation to get SNR in Spectral domain
if(exist('SNR','var') && SNR > 0)
    SNR_spectrum = SNR / sqrt(numel(t));
    std_FID = S_0/(2*SNR_spectrum);
else
    std_FID = 0;
end

% Create gaußian noise with std std_FID and mean 0
Noise = std_FID*complex(randn([1 vecSize]), randn([1 vecSize])); % Is real and imaginary noise uncorrelated?




%% Create FIDs

FID = [t; S_0 * exp(-(omega - LarmorFreq*2*pi) * 1i*t).*exp(-t/T2)*exp(1i*deg2rad(phase0)) + Noise];

%% Smooth FID

% if(exist('SmoothFIDSpan','var') && ne(SmoothFIDSpan,1))
%    t = AcqDelay:SmoothFIDSpan*dwelltime:t_end;
%    FID_smooth = transpose(smooth(FID(2,:),SmoothFIDSpan));
%    FID = [t; FID_smooth(1:SmoothFIDSpan:end)]; 
% end

if(exist('SmoothFIDSpan','var') && ne(SmoothFIDSpan,1))
   FID_smooth = transpose(smooth(FID(2,:),SmoothFIDSpan));
   FID = [t; FID_smooth]; 
   %FID = FID(:,1:end-SmoothFIDSpan);
end





%% Compute PPM-Scale

Chemshift = compute_chemshift_vector_1_2(LarmorFreq,dwelltime,numel(FID(2,:)));


%% FFT
Spectrum = [Chemshift; fftshift(fft(FID(2,:))) / sqrt(numel(t))];


end

