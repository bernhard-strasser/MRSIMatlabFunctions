function Noise_std = noise_computation_0_2(Signal,NoiseRegion,PointsPerSubregion,MeanThreshold)
%
% Compute the noise of a signal (std of the signal). It divides the signal in several subregions. If there is a higher std in a certain subregion than in the others it will be excluded.
% If the standard deviations are varying very severe (if the standard deviation of all the standard deviations is high) --> abort program with error
%
%

%% 0. Declarations, Preparations, Definitions

% Declarations


% Definitions
    

% 0.3 Preparations

% 0.3.1 Defining Subregions
TotalPeakregionPoints = NoiseRegion(2) - NoiseRegion(1) + 1;
TotalSubregions = floor(TotalPeakregionPoints/PointsPerSubregion);
ExtraPoints = mod(TotalPeakregionPoints,TotalSubregions);         % eg: have 403 TotalPeakregionPoints, 20 PointsPerSubregion --> 20 Subregions, and 3 extra point. Distribute these extrapoints to all SubRegions.


% 0.3.2 Compute a polyfit for real and imaginary part. Subtract this polyfit
Polyfit_Signal_Real = polyfit(NoiseRegion(1):NoiseRegion(2),real(Signal(NoiseRegion(1):NoiseRegion(2))),1);
Polyfit_Signal_Imag = polyfit(NoiseRegion(1):NoiseRegion(2),imag(Signal(NoiseRegion(1):NoiseRegion(2))),1);

Signal_SubBas = Signal;
Signal_SubBas(NoiseRegion(1):NoiseRegion(2)) = Signal(NoiseRegion(1):NoiseRegion(2)) - Polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)) - 1i*Polyval(Polyfit_Signal_Imag,NoiseRegion(1):NoiseRegion(2));


% 0.3.3 Initializations 
Mean_Subregions = zeros([1 TotalSubregions]);



%% 1. Compute the mean value of each subregion
for Subregion = 1:TotalSubregions-1
    Mean_Subregion(Subregion) = mean(Signal_SubBas((Subregion-1)*PointsPerSubregion + 1 : Subregion*PointsPerSubregion));
end
Mean_Subregion(TotalSubregion) = mean(Signal_SubBas((TotalSubregion-1)*PointsPerSubregion + 1 : NoiseRegion(2)));

Logical_MeanThreshold = Mean_Subregion > MeanThreshold;







