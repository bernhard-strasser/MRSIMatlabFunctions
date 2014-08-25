function [FID,Spectrum] = Simulate_FID_Spectra_1_0(omega,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,SmoothFIDSpan)
% Simulates FIDs. Input: Angular Frequency [rad/s], Zero-Order-Phase [Degree], AcquisitionDelay[s], T2[s], Magnitude of Signal at t = 0, SNR, dwelltime[s], vecSize, Span to Smooth FID

%% Comments, Definitions, Initializations
% Define Variables

water_frequency = 297*10^6; % @ 7T


% Define time
t_end = AcqDelay + dwelltime*(vecSize - 1);
t=AcqDelay:dwelltime:t_end;

% Define Standard Deviation to get SNR
std_FID = S_0/(2*SNR);

% Create gaußian noise with std std_FID and mean 0
Noise = std_FID*randn([1 vecSize]);




%% Create FIDs

FID = [t; S_0 * exp(omega * 1i*t).*exp(-t/T2)*exp(1i*degtorad(phase0)) + Noise + 1i*Noise]; % Is real and imaginary noise uncorrelated?

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

Chemshift = compute_chemshift_vector_1_1(water_frequency,dwelltime,numel(FID(2,:)));


%% FFT
Spectrum = [Chemshift; fftshift(fft(FID(2,:)))];


end

