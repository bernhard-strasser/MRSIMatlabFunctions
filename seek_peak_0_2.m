function [peak_found,SNR,Peak_Point,LeftBasis_Point,RightBasis_Point,LeftBasis_Signal,RightBasis_Signal,Basis_Signal_Mean,LeftBasisRegion,RightBasisRegion] = ...
          seek_peak_0_2(Signal,SeekPeakRegion,Noise,seek_criteria,seek_criteria_values)
% 
% Seeks peak within a given region with criteria:
%  - 'SNR': checks if the assumed peak (max of region) is SNR_Criterion_value times higher than the border of the region; this value must be given to the function
%  - 'derivative':  checks if the derivative of the spectrum is derivative_Criterion_value times higher than at the border; 
%                   This value and a averaging value must be given (compute the derivative between neighbouring points, or farther points) NOT IMPLEMENTED YET

%% 0. Preparations

% 0.1 Search which criteria should be checked

SNR_crit_index = strcmpi(seek_criteria,'SNR_Criterion');
SNR_crit_flag = logical(sum(SNR_crit_index));

%SNR_crit_index = strfnd(seek_criteria,'SNR_Criterion');                         % strfind gives an array number for each string, where 'SNR_Criterion' occurs. If it doesn't occur it gives an empty array [].

%Signal_SeekPeakRegion = Signal(min(SeekPeakRegion):max(SeekPeakRegion));


%% 1. SNR Criterion

if(SNR_crit_flag)
    
    % 1.0 Get the Criterion Values
    SNR_threshold = seek_criteria_values{SNR_crit_index}(1);                    % SNR-threshold for considering some maximum as peak or not
    BasisRegion_Width = seek_criteria_values{SNR_crit_index}(2);                % Look in certain region for minimum (considered as Basis-point): This value gives the extent of the region
    LeftBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}(3);     % Look for Left Basis-point this value points away from the found peak 
    RightBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}(4);    % Same for right Basis-point
    
    
    % 1.1 Find Maximum of Signal
    Peak_Point = find(Signal(min(SeekPeakRegion):max(SeekPeakRegion)) == max(Signal(min(SeekPeakRegion):max(SeekPeakRegion))));
    Peak_Point = Peak_Point + min(SeekPeakRegion) - 1;
    
    
    % 1.2 Compute BasisRegions (regions where to search for minima)
    LeftBasisRegion = [Peak_Point - LeftBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point - LeftBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    RightBasisRegion = [Peak_Point + RightBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point + RightBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    
    
%     % 1.3 Find minimum in left and right BasisRegions
%     LeftBasis_Point = find(Signal(LeftBasisRegion(1):LeftBasisRegion(2)) == min(Signal(LeftBasisRegion(1):LeftBasisRegion(2))));
%     RightBasis_Point = find(Signal(RightBasisRegion(1):RightBasisRegion(2)) == min(Signal(RightBasisRegion(1):RightBasisRegion(2)))); 
%     % Because we only sought in the small RightBasisRegion
%     LeftBasis_Point = LeftBasis_Point + LeftBasisRegion(1) - 1;                                 
%     RightBasis_Point = RightBasis_Point + RightBasisRegion(1) - 1;

    LeftBasis_Point = round(mean([LeftBasisRegion(1),LeftBasisRegion(2)]));         % The LeftBasis_Signal is most probably nowhere in the spectrum, so just take the middle of the Region where to plot the LeftBasis_Signal
    RightBasis_Point = round(mean([RightBasisRegion(1),RightBasisRegion(2)]));
    
    LeftBasis_Signal = mean(Signal(LeftBasisRegion(1):LeftBasisRegion(2)));
    RightBasis_Signal = mean(Signal(RightBasisRegion(1):RightBasisRegion(2)));
    
    
    % 1.4 Compute Mean Basis Signal
    Basis_Signal_Mean = mean([LeftBasis_Signal,RightBasis_Signal]);
    
    
    % 1.5 Compute Signal, Compute SNR
    
    Signal = Signal(Peak_Point) - Basis_Signal_Mean;
    SNR = Signal / (2*Noise);
    
    
    % 1.6 Set SNR_criterion to 0 if no peak found, 1 if found; If no peak found, set SNR to NaN
    SNR_crit = logical(SNR > SNR_threshold);
    
    if(~SNR_crit)
        SNR = NaN;
    end
    
    
    
else
	SNR_crit = true; 			% If it gets not checked it should not have any influence for peak_found
end


%% 2. More criteria needed ???



%% 3. Compute the logical peak_found


peak_found = SNR_crit;
