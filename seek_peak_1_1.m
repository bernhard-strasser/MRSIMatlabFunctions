function [NoPeakFound,Peak_Point,Peak_Signal,Basis_Point,Basis_Signal,Basis_Signal_Mean,LeftBasisRegion,RightBasisRegion,Polyfit_Signal,Signal_SubBas] = ...
          seek_peak_1_1(Signal,SeekPeakRegion,Use_RealAbs,PolyfitOrder,Phantom_flag,seek_criteria,seek_criteria_values)
% 
% Seeks peak within a given region with criteria:
%  - 'SNR': checks if the assumed peak (max of region) is SNR_Criterion_value times higher than the border of the region; this value must be given to the function
%  - 

%% 0. Declarations, Preparations, Definitions

% Declarations

% Definitions


% Preparations


% 0.1 Search which criteria should be checked

SNR_crit_index = strcmpi(seek_criteria,'SNR_Criterion');
SNR_crit_flag = logical(sum(SNR_crit_index));

mean_crit_index = strcmpi(seek_criteria,'Mean_Criterion');
mean_crit_flag = logical(sum(mean_crit_index));

findpeaks_crit_index = strcmpi(seek_criteria,'findpeaks');
findpeaks_crit_flag = logical(sum(findpeaks_crit_index));

SubBas_crit_index = strcmpi(seek_criteria,'SubBas');
SubBas_crit_flag = logical(sum(SubBas_crit_index));


% 0.2
Signal_RealAbs = eval([ Use_RealAbs '(Signal);' ]);


% 0.3 Subtract baseline of real and imaginary spectra

% 0.3.1 Fit Signal
if(exist('PolyfitOrder','var') && PolyfitOrder > 0)
    Polyfit_Signal = polyfit(SeekPeakRegion(1):SeekPeakRegion(2), Signal_RealAbs(SeekPeakRegion(1):SeekPeakRegion(2)),PolyfitOrder);                       % Fit baseline
    %Signal_SubBas = Signal_RealAbs;
    %Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)) = Signal_RealAbs(SeekPeakRegion(1):SeekPeakRegion(2)) - polyval(Polyfit_Signal,SeekPeakRegion(1):SeekPeakRegion(2));
    Signal_SubBas = Signal_RealAbs - polyval(Polyfit_Signal,1:size(Signal,2));   
else
    Signal_SubBas = Signal_RealAbs;
    Polyfit_Signal = 0;
end

% Signal_Time = ifft(ifftshift(Signal));
% Signal_Time_Apod = Signal_Time .* exp(-ApodValue*((1:numel(Signal_Time))-1)/numel(Signal_Time));
% Signal_Baseline = fftshift(fft(Signal_Time_Apod));
% Signal_SubBas = eval([ Use_RealAbs '(Signal - Signal_Baseline);' ]);



% figure
% plot(abs(Signal(SeekPeakRegion(1):SeekPeakRegion(2))))
% hold on
% plot(polyval(Polyfit_Signal_Abs,SeekPeakRegion(1):SeekPeakRegion(2)),'r')
% hold off
% 
% figure
% plot(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)))
% waitforbuttonpress
% waitforbuttonpress



%% 1. SNR Criterion

if(SNR_crit_flag)
    
    % 1.0 Get the Criterion Values
    Noise = seek_criteria_values{SNR_crit_index}{1};                            % SNR-threshold for considering some maximum as peak or not
    SNR_threshold = seek_criteria_values{SNR_crit_index}{2};                    % SNR-threshold for considering some maximum as peak or not
    BasisRegion_Width = seek_criteria_values{SNR_crit_index}{3};                % Look in certain region for minimum (considered as Basis-point): This value gives the extent of the region
    LeftBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}{4};     % Look for Left Basis-point this value points away from the found peak 
    RightBasisRegion_PeakDistance = seek_criteria_values{SNR_crit_index}{5};    % Same for right Basis-point
    
    
    % 1.1 Find Maximum of Signal
    %Peak_Point = eval([ 'find(' Use_RealAbs '(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2))) == max(' Use_RealAbs '(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)))));' ]);
    Peak_Point = find(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)) == max(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2))));
    Peak_Point = Peak_Point + SeekPeakRegion(1) - 1;

    
    
    % 1.2 Compute BasisRegions
    LeftBasisRegion = [Peak_Point - LeftBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point - LeftBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    RightBasisRegion = [Peak_Point + RightBasisRegion_PeakDistance - round(BasisRegion_Width/2), Peak_Point + RightBasisRegion_PeakDistance + round(BasisRegion_Width/2)];
    
    
%     % 1.3 Find minimum in left and right BasisRegions
%     LeftBasis_Point = find(Signal_SubBas(LeftBasisRegion(1):LeftBasisRegion(2)) == min(Signal_SubBas(LeftBasisRegion(1):LeftBasisRegion(2))));
%     RightBasis_Point = find(Signal_SubBas(RightBasisRegion(1):RightBasisRegion(2)) == min(Signal_SubBas(RightBasisRegion(1):RightBasisRegion(2)))); 
%     % Because we only sought in the small RightBasisRegion
%     LeftBasis_Point = LeftBasis_Point + LeftBasisRegion(1) - 1;                                 
%     RightBasis_Point = RightBasis_Point + RightBasisRegion(1) - 1;

    LeftBasis_Point = round(mean([LeftBasisRegion(1),LeftBasisRegion(2)]));         % The LeftBasis_Signal is most probably nowhere in the spectrum, so just take the middle of the Region where to plot the LeftBasis_Signal
    RightBasis_Point = round(mean([RightBasisRegion(1),RightBasisRegion(2)]));
    
    
%     % 1.4 Fit the Signal in the region from the one basis Region to the other
%     if(exist('PolyfitOrder','var') && PolyfitOrder > 0)
%         Polyfit_Signal = polyfit(LeftBasisRegion(1):RightBasisRegion(2), Signal_RealAbs(LeftBasisRegion(1):RightBasisRegion(2)),PolyfitOrder);                       % Fit baseline
%         Signal_SubBas = Signal_RealAbs;
%         Signal_SubBas(LeftBasisRegion(1):RightBasisRegion(2)) = Signal_RealAbs(LeftBasisRegion(1):RightBasisRegion(2)) - polyval(Polyfit_Signal,LeftBasisRegion(1):RightBasisRegion(2));
%     else
%         Signal_SubBas = Signal_RealAbs;
%         Polyfit_Signal = 0;
%     end       

    
    
    
%     % 1.4 Compute Signals
%     Peak_Signal = Signal_RealAbs(Peak_Point);
%     LeftBasis_Signal = mean(Signal_RealAbs(LeftBasisRegion(1):LeftBasisRegion(2)));
%     RightBasis_Signal = mean(Signal_RealAbs(RightBasisRegion(1):RightBasisRegion(2)));
%     Basis_Signal_Mean = mean([LeftBasis_Signal,RightBasis_Signal]);

    % 1.4 Compute Signals
    Peak_Signal = Signal_SubBas(Peak_Point);
    LeftBasis_Signal = mean(Signal_SubBas(LeftBasisRegion(1):LeftBasisRegion(2)));
    RightBasis_Signal = mean(Signal_SubBas(RightBasisRegion(1):RightBasisRegion(2)));
    if(Phantom_flag)
        %Valley_Point = find(Signal_SubBas(LeftBasisRegion(1):RightBasisRegion(2)) == min(Signal_SubBas(LeftBasisRegion(1):RightBasisRegion(2))));
        Basis_Signal_Mean = min(Signal_SubBas(LeftBasisRegion(1):RightBasisRegion(2)));
    else
        Basis_Signal_Mean = mean([LeftBasis_Signal,RightBasis_Signal]);
    end  
    
    % 1.5 Compute PeakHeight, Compute SNR    
    Peak_Height = Peak_Signal - Basis_Signal_Mean;
  
    
    
    
    SNR = Peak_Height / (2*Noise);
    
    
    % 1.6 Set SNR_criterion to 0 if no peak found, 1 if found; If no peak found, set SNR to NaN
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




%% 4. SubBas Criterion

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



%% DEBUG MODE: Plot


% figure
% plot(Signal_RealAbs(SeekPeakRegion(1):SeekPeakRegion(2)))
% if(PolyfitOrder > 0)
%     hold on  
%     %plot(real(Signal_Baseline(SeekPeakRegion(1):SeekPeakRegion(2))),'r')
%     plot(polyval(Polyfit_Signal,SeekPeakRegion(1):SeekPeakRegion(2)),'r')
%     hold off
% end
% 
% figure
% plot(Signal_SubBas(SeekPeakRegion(1):SeekPeakRegion(2)))
% figure
% plot(derivat,'r')
% hold on
% plot(derivat2,'g')
% hold off
% waitforbuttonpress

