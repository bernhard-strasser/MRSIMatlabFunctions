function [DensityEfficiencies,kSpaceDensities,kxPartitions_Ctr] = CalcTrajDensEff_2D(kSpaceTraj,BinSize)
%
% CalcEPSIEffWithRampSamp Calculate efficiency of EPSI trajectory
%
% This function was written by Bernhard Strasser, December 2017.
%
% The function takes as input the gradient values of the EPSI trajectory and the spectral bandwidth, and calculates the efficiency when doing
% ramp sampling due to the density difference to certain pre-defined target densities. Those target densities are defined now as
% Uniform, Hamming, and Gauss filter.
%
% [DensityEfficiencies,kSpaceDensities,kxPartitions_Ctr] = CalcEPSIEffWithRampSamp(GradX)
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






%% Calculate the density

% % Aproach with kSpaceDensity = 1/velocity:
% kSpaceDensity_1 = (1./GradX(2,1:end/2));

% VORONOI aproach: Calculating how many points are in equi-sized bins in (kx,ky). For this aproach, use a very small GRAD_RASTER_TIME (e.g. 10^-8).


kx_half = kSpaceTraj(1,:);
kxPartitions = min(kx_half):(max(kx_half)-min(kx_half))/BinSize:max(kx_half);
ky_half = kSpaceTraj(2,:);
kyPartitions = min(ky_half):(max(ky_half)-min(ky_half))/BinSize:max(ky_half);


NumSamplInBin = zeros([numel(kxPartitions)-1 numel(kyPartitions)-1]);
for jj=1:numel(kyPartitions)-1
    for ii=1:numel(kxPartitions)-1
        NumSamplInBin(ii,jj) = sum(kx_half >= kxPartitions(ii) & kx_half < kxPartitions(ii+1) & ky_half >= kyPartitions(jj) & ky_half < kyPartitions(jj+1));
    end
end
kSpaceDensities.Traj = NumSamplInBin/numel(kx_half)/numel(ky_half);

% Rescale the density, so that integral(kSpaceDensities.Traj,dk) = 1

% kXPartitions has length 1 larger than kSpaceDensities.Traj. Thus trapz( kXPartitions,kSpaceDensities.Traj) doesn't work. But we can define our partitons in the center of the bins
% But probably easier is the other way, to simply trapz(kSpaceDensities.Traj)*(kXPartitions(2)-kXPartitions(1)). Therefore this is commented out.
% kxPartitions_Ctr = kxPartitions(1:end-1) + kxPartSize/2;                   % Always the centers of the bins, these are one less

% Rescale so that integral(kSpaceDensities.Traj) = 1;
kxPartSize = kxPartitions(2)-kxPartitions(1);
kyPartSize = kyPartitions(2)-kyPartitions(1);

kSpaceDensities.Traj = kSpaceDensities.Traj/(trapz(trapz(kSpaceDensities.Traj))*kxPartSize*kyPartSize);

kSpaceDensities.Traj(kSpaceDensities.Traj == 0) = NaN;




%% Calculate the target densities

% Define and rescale hamming
kSpaceDensities.Hamm = HammingFilter(ones(size(kSpaceDensities.Traj,1)),[1 2],100,'OuterProduct',1);
kSpaceDensities.Hamm = EllipticalFilter(kSpaceDensities.Hamm,[1 2],[1 1 1 size(kSpaceDensities.Traj)/2]);
kSpaceDensities.Hamm = kSpaceDensities.Hamm/(trapz(trapz(kSpaceDensities.Hamm))*kxPartSize*kyPartSize);

% Define and rescale uniform density
kSpaceDensities.Uni =  ones(size(kSpaceDensities.Traj));
kSpaceDensities.Uni = EllipticalFilter(kSpaceDensities.Uni,[1 2],[1 1 1 size(kSpaceDensities.Traj)/2]);
kSpaceDensities.Uni = kSpaceDensities.Uni/(trapz(trapz(kSpaceDensities.Uni))*kxPartSize*kyPartSize);

% % Define and rescale Gauss
% kSpaceDensities.Gauss = 1:numel(kSpaceDensities.Traj); kSpaceDensities.Gauss = kSpaceDensities.Gauss - numel(kSpaceDensities.Gauss)/2 - 1; kSpaceDensities.Gauss = exp(-kSpaceDensities.Gauss.^2/80);
% kSpaceDensities.Gauss = kSpaceDensities.Gauss/(trapz(trapz(kSpaceDensities.Gauss))*kxPartSize*kyPartSize);

% kSpaceDensities.Traj = kSpaceDensities.Uni;


%% Calculate the efficiency based on used density and target density

kSpaceRatio = kSpaceDensities.Hamm.^2 ./ kSpaceDensities.Traj; kSpaceRatio(isnan(kSpaceRatio) | isinf(kSpaceRatio)) = 0;
DensityEfficiencies.VsHamm = trapz(trapz(kSpaceRatio))*kxPartSize*kyPartSize;
DensityEfficiencies.VsHamm = sqrt(1/DensityEfficiencies.VsHamm);

kSpaceRatio = kSpaceDensities.Uni.^2 ./ kSpaceDensities.Traj; kSpaceRatio(isnan(kSpaceRatio) | isinf(kSpaceRatio)) = 0;
DensityEfficiencies.VsUni = trapz(trapz(kSpaceRatio))*kxPartSize*kyPartSize;
DensityEfficiencies.VsUni = sqrt(1/DensityEfficiencies.VsUni);

% DensityEfficiencies.VsGauss = trapz(trapz(kSpaceDensities.Gauss.^2 ./ (kSpaceDensities.Traj)))*kxPartSize*kyPartSize;
% DensityEfficiencies.VsGauss = sqrt(1/DensityEfficiencies.VsGauss);





