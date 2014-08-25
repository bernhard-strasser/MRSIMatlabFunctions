function [NoPeakFound,Peak_Point,Peak_Signal,Basis_Point,Basis_Signal,Basis_Signal_Mean,LeftBasisRegion,RightBasisRegion,Polyfit_Signal_Real,Polyfit_Signal_Imag] = ...
          seek_peak_0_3(Signal,SeekPeakRegion,seek_criteria,seek_criteria_values)
% 
% Seeks peak within a given region with criteria:
%  - 'SNR': checks if the assumed peak (max of region) is SNR_Criterion_value times higher than the border of the region; this value must be given to the function
%  - 

%% 0. Preparations

% 0.1 Search which criteria should be checked

SNR_crit_index = strcmpi(seek_criteria,'SNR_Criterion');
SNR_crit_flag = logical(sum(SNR_crit_index));

mean_crit_index = strcmpi(seek_criteria,'Mean_Criterion');
mean_crit_flag = logical(sum(mean_crit_index));

findpeaks_crit_index = strcmpi(seek_criteria,'findpeaks');
findpeaks_crit_flag = logical(sum(findpeaks_crit_index));



%% 1. SNR Criterion

if(SNR_crit_flag)
    
    % 1.0 Get the Criterion Values
    Noise = seek_criteria_values{SNR_crit_index}{1};                            % SNR-threshold for considering some maximum as peak or not
    SNR_threshold = seek_criteria_values{SNR_crit_index}{2};                    % SNR-threshold for considering some maximum as peak or not
    BasisRegion_Width = seek_criteria_values{SNR_crit_index}{3};                % Look in certain region for minimum (considered as Basis-point): This value gives the extent of the region
    LeftBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}{4};     % Look for Left Basis-point this value points away from the found peak 
    RightBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}{5};    % Same for right Basis-point
    Use_RealAbs = seek_criteria_values{SNR_crit_index}{6};
    
    
    % 1.1 Find Maximum of Magnitude Signal
    Peak_Point = find(abs(Signal(SeekPeakRegion(1):SeekPeakRegion(2))) == max(abs(Signal(SeekPeakRegion(1):SeekPeakRegion(2)))));
    Peak_Point = Peak_Point + SeekPeakRegion(1) - 1;
    
    
    % 1.2 Compute BasisRegions
    LeftBasisRegion = [Peak_Point - LeftBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point - LeftBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    RightBasisRegion = [Peak_Point + RightBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point + RightBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    
    
%     % 1.3 Find minimum in left and right BasisRegions
%     LeftBasis_Point = find(Signal_RealAbs(LeftBasisRegion(1):LeftBasisRegion(2)) == min(Signal_RealAbs(LeftBasisRegion(1):LeftBasisRegion(2))));
%     RightBasis_Point = find(Signal_RealAbs(RightBasisRegion(1):RightBasisRegion(2)) == min(Signal_RealAbs(RightBasisRegion(1):RightBasisRegion(2)))); 
%     % Because we only sought in the small RightBasisRegion
%     LeftBasis_Point = LeftBasis_Point + LeftBasisRegion(1) - 1;                                 
%     RightBasis_Point = RightBasis_Point + RightBasisRegion(1) - 1;

    LeftBasis_Point = round(mean([LeftBasisRegion(1),LeftBasisRegion(2)]));         % The LeftBasis_Signal is most probably nowhere in the spectrum, so just take the middle of the Region where to plot the LeftBasis_Signal
    RightBasis_Point = round(mean([RightBasisRegion(1),RightBasisRegion(2)]));
    
    
    
    % 1.4 Subtract baseline of real and imaginary spectra

    % 1.4.1 Compute Region for which Baseline gets computed
    Baseline_Point = [min([LeftBasisRegion(1),SeekPeakRegion(1)]),max([RightBasisRegion(2),SeekPeakRegion(2)])];

    % 1.4.2 Fit Real Part
    Polyfit_Signal_Real = polyfit(Baseline_Point(1):Baseline_Point(2),real(Signal(Baseline_Point(1):Baseline_Point(2))),1);                       % Fit baseline
    Signal_Real_SubBas = real(Signal(Baseline_Point(1):Baseline_Point(2))) - polyval(Polyfit_Signal_Real,Baseline_Point(1):Baseline_Point(2));      % Subtract that baseline

    % 1.4.3 Fit Imaginary Part
    Polyfit_Signal_Imag = polyfit(Baseline_Point(1):Baseline_Point(2),imag(Signal(Baseline_Point(1):Baseline_Point(2))),1);                       % Fit baseline
    Signal_Imag_SubBas = imag(Signal(Baseline_Point(1):Baseline_Point(2))) - polyval(Polyfit_Signal_Imag,Baseline_Point(1):Baseline_Point(2));      % Subtract that baseline

    Signal_SubBas = Signal;
    Signal_SubBas(Baseline_Point(1):Baseline_Point(2)) = Signal_Real_SubBas + 1i*Signal_Imag_SubBas;
    
    
    
    % 1.5 Compute Signals
    Signal_RealAbs = eval([Use_RealAbs '(Signal_SubBas);']);
    Peak_Signal = Signal_RealAbs(Peak_Point);
    LeftBasis_Signal = mean(Signal_RealAbs(LeftBasisRegion(1):LeftBasisRegion(2)));
    RightBasis_Signal = mean(Signal_RealAbs(RightBasisRegion(1):RightBasisRegion(2)));
    Basis_Signal_Mean = mean([LeftBasis_Signal,RightBasis_Signal]);
    
    
    % 1.6 Compute PeakHeight, Compute SNR
    Peak_Height = Signal_RealAbs(Peak_Point) - Basis_Signal_Mean;
    SNR = Peak_Height / (2*Noise);
    
    
    % 1.7 Set SNR_criterion to 0 if no peak found, 1 if found; If no peak found, set SNR to NaN
    SNR_crit = logical(SNR > SNR_threshold);
    
    Peak_Point = {Peak_Point};
    Peak_Signal = {Peak_Signal};
    Basis_Point = {[LeftBasis_Point,RightBasis_Point]};
    Basis_Signal = {[LeftBasis_Signal,RightBasis_Signal]};
    Basis_Signal_Mean = {Basis_Signal_Mean};
    LeftBasisRegion = {LeftBasisRegion};
    RightBasisRegion = {RightBasisRegion};       
        
    NoPeakFound = ~SNR_crit;
          
    
end


%% 2. Mean Criterion

if(mean_crit_flag)

    % 2.0 Get the criterion values
    PointsPerSubregion = seek_criteria_values{mean_crit_index}(1);
    
    
    % 2.1 Preparations
    
    TotalPeakregionPoints = SeekPeakRegion(2) - SeekPeakRegion(1) + 1;
    TotalSubregions = round(TotalPeakregionPoints/PointsPerSubregion);
    
    % 2.2 Compute a polyfit for real and imaginary part
    
    Polyfit_Signal_RealAbs = polyfit(SeekPeakRegion(1):SeekPeakRegion(2),real(Signal(SeekPeakRegion(1):SeekPeakRegion(2))),3);
    Polyfit_Signal_Imag = polyfit(SeekPeakRegion(1):SeekPeakRegion(2),imag(Signal(SeekPeakRegion(1):SeekPeakRegion(2))),3);
    
    Signal_WithoutBaseline = Signal - Polyval(Polyfit_Signal_RealAbs,SeekPeakRegion(1):SeekPeakRegion(2)) - 1i*Polyval(Polyfit_Signal_Imag,SeekPeakRegion(1):SeekPeakRegion(2));
    
    
    % 2.3 - 2.8 while loop: search for peaks as long as peaks are found
    search_on = true;
    while(search_on)
       
        % 2.3 Find maximum in magnitude spectrum
        Peak_Point = find(abs(Signal_WithoutBaseline(min(SeekPeakRegion):max(SeekPeakRegion))) == max(abs(Signal_WithoutBaseline(min(SeekPeakRegion):max(SeekPeakRegion)))));
        Peak_Point = Peak_Point + min(SeekPeakRegion) - 1;
        
        % 2.4 Divide the SeekPeakRegion in Subregions with PointsPerSubregion points

        Subregions 
        
    end


else
    mean_crit = true;
end


%% 3. findpeaks Function

if(findpeaks_crit_flag)
    
    % 3.0 Get the Criterion Values
    SlopeThreshold = seek_criteria_values{finpeaks_crit_index}{1};                            % Threshold for first derivative: If 
    AmpThreshold = seek_criteria_values{finpeaks_crit_index}{2};                    % SNR-threshold for considering some maximum as peak or not
    smoothwidth = seek_criteria_values{finpeaks_crit_index}{3};                % Look in certain region for minimum (considered as Basis-point): This value gives the extent of the region
    peakgroup = seek_criteria_values{finpeaks_crit_index}{4};     % Look for Left Basis-point this value points away from the found peak 
    smoothtype = seek_criteria_values{finpeaks_crit_index}{5};    % Same for right Basis-point


    Peaks = findpeaks(abs(Signal_SubBas),PeakRegion(1):PeakRegion(2),SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
    

    
    
    Peak_Point = Peaks + Noise_Peak_SP - 1;
    Noise_PeakRegion_HiFi_SP = Noise_Peak_SP - round(Noise_Peak_Width/2);  
    Noise_PeakRegion_LoFi_SP = Noise_Peak_SP + round(Noise_Peak_Width/2); 

end

