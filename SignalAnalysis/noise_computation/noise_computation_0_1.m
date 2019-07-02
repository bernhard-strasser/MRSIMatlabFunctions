function Noise_std = noise_computation_0_1(Signal,NoiseRegion,TotalSubregions)
%
% Compute the noise of a signal (std of the signal). It divides the signal in several subregions. If there is a higher std in a certain subregion than in the others it will be excluded.
% If the standard deviations are varying very severe (if the standard deviation of all the standard deviations is high) --> abort program with error
%
%

%% 0. Declarations, Preparations, Definitions

% Declarations


% Definitions

MinimumPointsPerSubregion = 50;
    

% Preparations

TotalPoints = NoiseRegion(2) - NoiseRegion(1) + 1;
if(TotalPoints < MinimumPointsPerSubregion)
    display([ char(10) 'Error: Not enough points within the Noise Region. Aborting program.'])
    return;
end
PointsPerSubregion = round(TotalPoints/TotalSubregions);    % The first 1:TotalSubregions-1 Subregions have this number of points, the last one has PointsPerSubregion + all the remaining points (remaining because of round)

while(PointsPerSubregion < MinimumPointsPerSubregion)
    display([char(10) 'Error: Too many subregions. Decreasing TotalSubregions.' ])
    TotalSubregions = TotalSubregions - 1;
    PointsPerSubregion = round(TotalPoints/TotalSubregions);
end


Subregions_Standarddevs = zeros([1 TotalSubregions]);



%% 1. Compute sum of signal minus linear fit for each subregion
% The sum of gaussian noise is approximately zero. If 
for k=1:TotalSubregions-1                                                                                            % The last SubRegion is treated different 
    Subregions_Standarddevs(1) = std(Signal(NoiseRegion(1) + (k-1)*PointsPerSubregion : NoiseRegion(1) + k*PointsPerSubregion));
end
Subregions_Standarddevs(TotalSubregions) = std(Signal(NoiseRegion(1) + (k+1)*PointsPerSubregion : NoiseRegion(2)));  % The last SubRegion has all the remain points inside


