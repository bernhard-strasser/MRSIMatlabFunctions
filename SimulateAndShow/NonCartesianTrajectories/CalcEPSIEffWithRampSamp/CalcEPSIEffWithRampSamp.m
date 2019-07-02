function [DensityEfficiencies,kSpaceDensities,kxPartitions_Ctr] = CalcEPSIEffWithRampSamp(GradX,SpectralBW)
%
% CalcEPSIEffWithRampSamp Calculate efficiency of EPSI trajectory
%
% This function was written by Bernhard Strasser, December 2017.
%
% The function takes as input the gradient values of the EPSI trajectory and the spectral bandwidth, and calculates the efficiency when doing
% ramp sampling due to the density difference to certain pre-defined target densities. Those target densities are defined now as
% Uniform, Hamming, and Gauss filter.
%
% [DensityEfficiencies,kSpaceDensities,kxPartitions_Ctr] = CalcEPSIEffWithRampSamp(GradX,SpectralBW)
%
% Input: 
% -         GradX                             ...    Matrix of size [2 N] storing the time values (GradX(1,:) in s (?), and the gradient values in
%                                                    x-direction (GradX(2,:)) in T/m.
% -         SpectralBW                        ...    Spectral Bandwidth
% Output:
% -         DensityEfficiencies               ...    The calculated efficiencies in comparison to the pre-defined densities.
% -         kSpaceDensities                   ...    The calculated k-space densities, useful for plotting.
% -         kxPartitions_Ctr                  ...    The center of the equi-distant bins ( = Partitions) in which we divide our kx-values.
%                                                    Useful for plotting the kSpaceDensities.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ???

% Further remarks: 

%% Definitions


% if(exist('Par','var') && ~isstruct(Par))
%     fprintf('\n\nWarning: The input variable ''Par'' needs to be a structure, see help CalculateEPSITrajectory. Assign standard values.\n')
%     clear Par;
% end
%        
%     
% if(~exist('Par','var') || ~isfield(Par,'Sys') || ~isstruct(Par.Sys))
%     Par = rmfield(Par,'Sys');
%     Par.Sys.Gamma = 42.5756*10^6;
%     Par.Sys.SlewRate_Max = 200;
%     Par.Sys.GRAD_RASTER_TIME = 1E-6;
% end
% if(~isfield(Par.Sys,'Gamma') || Par.Sys.Gamma == 0)
%     Par.Sys.Gamma = 42.5756*10^6;
% end
% if(~isfield(Par.Sys,'SlewRate_Max') || Par.Sys.SlewRate_Max == 0)
%     Par.Sys.SlewRate_Max = 200;
% end
% if(~isfield(Par.Sys,'GRAD_RASTER_TIME') || Par.Sys.GRAD_RASTER_TIME == 0)
%     Par.Sys.GRAD_RASTER_TIME = 1E-6;
% end
% 
% if(~isfield(Par,'MeasProt') || ~isstruct(Par.MeasProt))
%     Par = rmfield(Par,'MeasProt');
%     Par.MeasProt.RampSampling_Flag = true;    
%     Par.MeasProt.SpectrRange_Target = 2700;
%     Par.MeasProt.nTI = 1;
%     Par.MeasProt.Nx = 8;
%     Par.MeasProt.FOVx = 0.22;
%     Par.MeasProt.vecSizePerTI = 4;
% end
% if(~isfield(Par.MeasProt,'RampSampling_Flag'))
%     Par.MeasProt.RampSampling_Flag = true;
% end
% if(~isfield(Par.MeasProt,'SpectrRange_Target') || Par.MeasProt.SpectrRange_Target == 0)
%     Par.MeasProt.SpectrRange_Target = 2700;
% end
% if(~isfield(Par.MeasProt,'nTI') || Par.MeasProt.nTI == 0)
%     Par.MeasProt.nTI = 1;
% end
% if(~isfield(Par.MeasProt,'Nx') || Par.MeasProt.Nx == 0)
%     Par.MeasProt.Nx = 8;
% end
% if(~isfield(Par.MeasProt,'FOVx') || Par.MeasProt.FOVx == 0)
%     Par.MeasProt.FOVx = 0.22;
% end
% if(~isfield(Par.MeasProt,'vecSizePerTI') || Par.MeasProt.vecSizePerTI == 0)
%     Par.MeasProt.vecSizePerTI = 4;
% end






%% Calculate the density

% % Aproach with kSpaceDensity = 1/velocity:
% kSpaceDensity_1 = (1./GradX(2,1:end/2));

% VORONOI aproach: Calculating how many points are in equi-sized bins along kx. For this aproach, use a very small GRAD_RASTER_TIME (e.g. 10^-8).

kx = cumtrapz(GradX(1,:),GradX(2,:));

ZeroCrossing = find(abs(GradX(2,:)) < 10^-8);
kx = kx - kx(ZeroCrossing(2))/2;
kx_half = kx;
kxPartitions = min(kx_half):(max(kx_half)-min(kx_half))/64:max(kx_half);
NumSamplInBin = zeros([1 numel(kxPartitions)-1]);
for i=1:numel(kxPartitions)-1
    NumSamplInBin(i) = sum(kx_half >= kxPartitions(i) & kx_half < kxPartitions(i+1));
end
kSpaceDensities.EPSI = NumSamplInBin/numel(kx_half);

% Rescale the density, so that integral(kSpaceDensities.EPSI,dk) = 1/SpectralBW

% kXPartitions has length 1 larger than kSpaceDensities.EPSI. Thus trapz( kXPartitions,kSpaceDensities.EPSI) doesn't work. But we can define our partitons in the center of the bins
% But probably easier is the other way, to simply trapz(kSpaceDensities.EPSI)*(kXPartitions(2)-kXPartitions(1)). Therefore this is commented out.
kxPartitions_Ctr = kxPartitions(1:end-1) + (kxPartitions(2)-kxPartitions(1))/2;                   % Always the centers of the bins, these are one less

% Rescale so that integral(kSpaceDensities.EPSI) = 1;
kSpaceDensities.EPSI = kSpaceDensities.EPSI/(trapz(kSpaceDensities.EPSI)*(kxPartitions(2)-kxPartitions(1)));

% Rescale so that integral(kSpaceDensities.EPSI) = 1/SpectralBW
kSpaceDensities.EPSI = kSpaceDensities.EPSI/SpectralBW;







%% Calculate the target densities

% Define and rescale hamming
kSpaceDensities.Hamm = transpose(hamming(numel(kSpaceDensities.EPSI)));
kSpaceDensities.Hamm = kSpaceDensities.Hamm/(trapz(kSpaceDensities.Hamm)*(kxPartitions(2)-kxPartitions(1))*SpectralBW);

% Define and rescale uniform density
kSpaceDensities.Uni = ones(size(kSpaceDensities.EPSI));
kSpaceDensities.Uni = kSpaceDensities.Uni/(trapz(kSpaceDensities.Uni)*(kxPartitions(2)-kxPartitions(1))*SpectralBW);

% Define and rescale Gauss
kSpaceDensities.Gauss = 1:numel(kSpaceDensities.EPSI); kSpaceDensities.Gauss = kSpaceDensities.Gauss - numel(kSpaceDensities.Gauss)/2 - 1; kSpaceDensities.Gauss = exp(-kSpaceDensities.Gauss.^2/80);
kSpaceDensities.Gauss = kSpaceDensities.Gauss/(trapz(kSpaceDensities.Gauss)*(kxPartitions(2)-kxPartitions(1))*SpectralBW);

% kSpaceDensities.EPSI = kSpaceDensities.Uni;


%% Calculate the efficiency based on used density and target density

DensityEfficiencies.VsHamm = trapz(kSpaceDensities.Hamm.^2 ./ (kSpaceDensities.EPSI)) * (kxPartitions(2)-kxPartitions(1)) * SpectralBW;
DensityEfficiencies.VsHamm = sqrt(1/DensityEfficiencies.VsHamm);

DensityEfficiencies.VsUni = trapz(kSpaceDensities.Uni.^2 ./ (kSpaceDensities.EPSI)) * (kxPartitions(2)-kxPartitions(1)) * SpectralBW;
DensityEfficiencies.VsUni = sqrt(1/DensityEfficiencies.VsUni);

DensityEfficiencies.VsGauss = trapz(kSpaceDensities.Gauss.^2 ./ (kSpaceDensities.EPSI)) * (kxPartitions(2)-kxPartitions(1)) * SpectralBW;
DensityEfficiencies.VsGauss = sqrt(1/DensityEfficiencies.VsGauss);





