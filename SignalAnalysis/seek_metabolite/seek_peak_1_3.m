function [NoPeakFound,Peak_Point,Peak_Signal,Basis_Point,Basis_Signal,Basis_Signal_Mean,LeftBasisRegion,RightBasisRegion,Polyfit_Signal,Signal_SubBas] = ...
          seek_peak_1_3(Signal,SeekPeakRegion,Use_RealAbs,PolyfitOrder,Phantom_flag,seek_criteria,seek_criteria_values)
%
% seek_peak_x_y, find peaks within a spectrum with different criteria.
%
% This function was written by Bernhard Strasser, 2011 - 2012.
%
%
% The function searches for a peak and decides with different criteria if it is a peak or not. It then gives you where the peak and where the basis of the peak was found, and other stuff.
% The criteria are:
% - 'SNR_Criterion': checks if the assumed peak (max of region) is SNR_Criterion_value times higher than the border of the region; this value must be given to the function
% - No other working criteria yet, only ideas.
%
% [SNR_mat,Shift_mat] = compute_SNR_matrix_x_y(InData,OutInfo,LarmorFreq,dwelltime,zerofilling,NoiseRegion_CS,SeekMetabo,Phantom_flag)
%
%
% Input:
% Signal:               The signal in frequency domain.  
%
% SeekPeakRegion:       The region in points where the peak should be searched for in the signal, e.g. SeekPeakRegion = [30, 60]
% Use_RealAbs:          If the peak should be searched in real or absolute spectrum.
% PolyfitOrder:         The order of the polynom. This polynom defines the "Baseline" of the spectrum and gets subtracted. If PolyfitOrder = 0, nothing is done. 
%                       This definition of the Baseline lead to big troubles, avoid using it.
% Phantom_flag:         If the signal comes from a phantom data set. The "Basis Points" are computed differently in that case.
% seek_criteria:        Cell array of character arrays determining which peak-criteria should be used. Several criteria can be used independently for checking for the peak.
%                       Possible entries: 'SNR_Criterion', 'Mean_Criterion' (not working!), 'findpeaks' (not working!), 'SubBas' (not working!).
% seek_criteria_values: Cell array of cell arrays. For each criterion you want to use you have to pass over a cell array of values, bundled together in a cell array. These values are
%                       used to control the peak-finding process.
%                       In case of the only working criterion, the 'SNR_Criterion', these values are:
%                           * seek_criteria_values{1}: Noise:  2*standard deviation of the noise.
%                           * seek_criteria_values{2}: SNR_threshold: SNR-threshold for considering some maximum as peak or not
%                           * seek_criteria_values{3}: SimilarBasisSignals_Ratio: Used for the SimilarBasisSignals-SubCriterion. See Below.
%                           * seek_criteria_values{4}: BasisRegion_Width: This value gives the extent of the regions where to search for the Basis points
%                           * seek_criteria_values{5}: LeftBasisRegion_PeakDistance: Distance of the left basis region to the found (or assumed) peak
%                           * seek_criteria_values{6}: RightBasisRegion_PeakDistance: Same for right Basis-point
%
%
% Output
% NoPeakFound:          If true, within the SeekPeakRegion nothing was found that fulfils the criteria for countin as a peak.
% Peak_Point:           The spectral point where the peak was found (or searched if no peak was found).
% Basis_Point:          The points defined as the basis of the peak to its left and right. The difference of the peak minus these points is the peak height.
% Basis_Signal:         The Signal of the spectrum at the Basis_Points
% Basis_Signal_Mean:    The mean of the right and left Basis_Signals
% LeftBasisRegion:      The region where the left Basis_Point was searched
% RightBasisRegion:     Like Left.
% Polyfit_Signal:       The polynom-spectrum.
% Signal_SubBas:        Signal - Polyfit_Signal
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks:




      % 
% Seeks peak within a given region with criteria:
%  - 'SNR_Criterion': checks if the assumed peak (max of region) is SNR_Criterion_value times higher than the border of the region; this value must be given to the function
%  - 



%% 0. Declarations, Preparations, Definitions

% Declarations

% Definitions


% Preparations



% 0.1 Search which criteria should be checked and in which cell-array-index of seek-criteria this can be found. Easier would be to use a structure.
% Several criteria can be checked for the same criteria, thus this is so complicated.

SNR_crit_index = strcmpi(seek_criteria,'SNR_Criterion');
SNR_crit_flag = logical(sum(SNR_crit_index));

mean_crit_index = strcmpi(seek_criteria,'Mean_Criterion');
mean_crit_flag = logical(sum(mean_crit_index));

findpeaks_crit_index = strcmpi(seek_criteria,'findpeaks');
findpeaks_crit_flag = logical(sum(findpeaks_crit_index));

SubBas_crit_index = strcmpi(seek_criteria,'SubBas');
SubBas_crit_flag = logical(sum(SubBas_crit_index));



% 0.3
Signal_RealAbs = eval([ Use_RealAbs '(Signal);' ]);


% 0.4 Subtract baseline of real and imaginary spectra

% 0.4.1 Fit Signal
if(exist('PolyfitOrder','var') && PolyfitOrder > 0)
    Polyfit_Signal = polyfit(SeekPeakRegion(1):SeekPeakRegion(2), Signal_RealAbs(SeekPeakRegion(1):SeekPeakRegion(2)),PolyfitOrder);                       % Fit baseline
    Signal_SubBas = Signal_RealAbs - polyval(Polyfit_Signal,1:size(Signal,2));   
else
    Signal_SubBas = Signal_RealAbs;
    Polyfit_Signal = 0;
end

% Instead of PolyFit, use apodized signal as Baseline
% Signal_Time = ifft(ifftshift(Signal));
% Signal_Time_Apod = Signal_Time .* exp(-ApodValue*((1:numel(Signal_Time))-1)/numel(Signal_Time));
% Signal_Baseline = fftshift(fft(Signal_Time_Apod));
% Signal_SubBas = eval([ Use_RealAbs '(Signal - Signal_Baseline);' ]);



%% 1. SNR Criterion

if(SNR_crit_flag)

    
    % 1.0 Get the Criterion Values
    % If input variables were omitted.
    if(numel(seek_criteria_values{SNR_crit_index}) < 4)
        BasisRegion_Width = 0;
        LeftBasisRegion_PeakDistance = 0;
        RightBasisRegion_PeakDistance = 0;
    else
        BasisRegion_Width = seek_criteria_values{SNR_crit_index}{4};                % Look in certain region for minimum (considered as Basis-point): This value gives the extent of the region
        LeftBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}{5};     % Look for Left Basis-point LeftBasisRegion_PeakDistance away from the found peak 
        RightBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}{6};    % Same for right Basis-point
    end

    
    Noise = seek_criteria_values{SNR_crit_index}{1};                                % 2*standard deviation of the noise.
    SNR_threshold = seek_criteria_values{SNR_crit_index}{2};                        % SNR-threshold for considering some maximum as peak or not
    SimilarBasisSignals_Ratio = seek_criteria_values{SNR_crit_index}{3};            % Used for the SimilarBasisSignals-SubCriterion. See Below.
    
    
    % 1.1 Find Maximum of Signal
    Peak_Point = find(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)) == max(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2))));
    Peak_Point = Peak_Point + SeekPeakRegion(1) - 1;

    
    % 1.2 Define BasisRegion, i.e. the region where to "search" for the "Basis", and the points "Left/RightBasis_Point" which define where th 
    LeftBasisRegion = [Peak_Point - LeftBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point - LeftBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    RightBasisRegion = [Peak_Point + RightBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point + RightBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    LeftBasis_Point = round(mean([LeftBasisRegion(1),LeftBasisRegion(2)]));         % The LeftBasis_Signal is most probably nowhere in the spectrum, so just take the middle of the Region where to plot the LeftBasis_Signal
    RightBasis_Point = round(mean([RightBasisRegion(1),RightBasisRegion(2)]));
 
    
    % 1.3 Compute Signals
    Peak_Signal = Signal_SubBas(Peak_Point);
    
    
    % 1.4 Compute the Basis Signal
    if(LeftBasisRegion_PeakDistance == 0 || RightBasisRegion_PeakDistance == 0)
        LeftBasis_Signal = 0;
        RightBasis_Signal = 0;
        Basis_Signal_Mean = 0;
        
    else
        LeftBasis_Signal = mean(Signal_SubBas(LeftBasisRegion(1):LeftBasisRegion(2)));
        RightBasis_Signal = mean(Signal_SubBas(RightBasisRegion(1):RightBasisRegion(2)));
        if(Phantom_flag)
            Basis_Signal_Mean = min(Signal_SubBas(LeftBasisRegion(1):RightBasisRegion(2)));
        else
            Basis_Signal_Mean = mean([LeftBasis_Signal,RightBasis_Signal]);
        end 
        
    end
    
    
    % 1.5 Compute PeakHeight, Compute SNR    
    Peak_Height = Peak_Signal - Basis_Signal_Mean;
    SNR = Peak_Height / (2*Noise);   
    
    
    % 1.6 Bais-Below-Peak Criterion: Both Basis Signals must be below the Peak.
    if(LeftBasis_Signal < Peak_Signal && RightBasis_Signal < Peak_Signal)
        BasisBelowPeak_crit = true;
    else
        BasisBelowPeak_crit = false;
    end

    
    % 1.7 SimilarBasisSignals Criterion: The Basis signals should be similar in height.
    if(abs(LeftBasis_Signal - RightBasis_Signal) < abs(Peak_Signal - max(LeftBasis_Signal,RightBasis_Signal))*SimilarBasisSignals_Ratio)     % 4/5; The higher the value the more tolerant this criterion is. 
        SimilarBasisSignals_crit = true;                                                                                                     % If e.g. this value is 1, the criterion means that the Basis_Signals must be  
    else                                                                                                                                     % closer than the Peak-signal and the higher Basis signal.
        SimilarBasisSignals_crit = false;
    end
    
    
    % 1.8 SNR_criterion
    SNR_crit = logical(SNR > SNR_threshold);
    
    
    % 1.9 Combine All Criteria
    NoPeakFound = ~(SNR_crit * BasisBelowPeak_crit * SimilarBasisSignals_crit);    


    % 1.10 Output Assignments
    Basis_Point = [LeftBasis_Point,RightBasis_Point];
    Basis_Signal = [LeftBasis_Signal,RightBasis_Signal];
   
  

    
end



%% 2. NOT WORKING: Mean Criterion

if(mean_crit_flag)

    % 2.0 Get the criterion values
    PointsPerSubregion = seek_criteria_values{mean_crit_index}{1};
    
    
    % 2.1 Preparations
    
    TotalPeakregionPoints = SeekPeakRegion(2) - SeekPeakRegion(1) + 1;
    TotalSubregions = round(TotalPeakregionPoints/PointsPerSubregion);
    
    
    % 2.2 Compute a polyfit for real and imaginary part and subtract that from the signal
    
    Polyfit_Signal_Real = polyfit(SeekPeakRegion(1):SeekPeakRegion(2),real(Signal(SeekPeakRegion(1):SeekPeakRegion(2))),1);
    Polyfit_Signal_Imag = polyfit(SeekPeakRegion(1):SeekPeakRegion(2),imag(Signal(SeekPeakRegion(1):SeekPeakRegion(2))),1);
    
    Signal_SubBas = Signal - Polyval(Polyfit_Signal_Real,SeekPeakRegion(1):SeekPeakRegion(2)) - 1i*Polyval(Polyfit_Signal_Imag,SeekPeakRegion(1):SeekPeakRegion(2));
    
    
    % 2.3 - 2.8 while loop: search for peaks as long as peaks are found
    search_on = true;
    while(search_on)
       
        % 2.3 Find maximum in spectrum
        Peak_Point = eval(['find(' Use_RealAbs '(Signal_SubBas(min(SeekPeakRegion):max(SeekPeakRegion))) == max(' Use_RealAbs '(Signal_SubBas(min(SeekPeakRegion):max(SeekPeakRegion)))));']);
        Peak_Point = Peak_Point + min(SeekPeakRegion) - 1;
        
        % 2.4 Divide the SeekPeakRegion in Subregions with PointsPerSubregion points
        Subregions 
        
    end


else
    mean_crit = true;
end


%% 3. NOT WORKING: findpeaks Function

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


%% 4. NOT WORKING: SubBas Criterion

if(SubBas_crit_flag)
    
    % 4.1 Find Maximum of Signal
    Peak_Point = find(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)) == max(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2))));
    Peak_Point = Peak_Point + SeekPeakRegion(1) - 1;
    
    % 4.2 Basis
    
    LeftBasisRegion = [Peak_Point - 15, Peak_Point - 5];
    RightBasisRegion = [Peak_Point + 5, Peak_Point + 15];
    LeftBasis_Point = Peak_Point - 10;
    RightBasis_Point = Peak_Point + 10;

    % 1.4 Compute Signals
    Peak_Signal = Signal_SubBas(Peak_Point);
    LeftBasis_Signal = 0;
    RightBasis_Signal = 0;
    Basis_Signal_Mean = 0;
    
    
    Peak_Point = {Peak_Point};
    Peak_Signal = {Peak_Signal};
    Basis_Point = {[LeftBasis_Point,RightBasis_Point]};
    Basis_Signal = {[LeftBasis_Signal,RightBasis_Signal]};
    Basis_Signal_Mean = {Basis_Signal_Mean};
    LeftBasisRegion = {LeftBasisRegion};
    RightBasisRegion = {RightBasisRegion}; 
      
    NoPeakFound = false;
    
end


