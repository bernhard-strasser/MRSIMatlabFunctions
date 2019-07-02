function [Noise_std,Signal_ApodSmooth,Signal_SubPeak,ExcludeSubregions_StartEndPoints] = noise_computation_new_0_5(Signal,NoiseRegion,ApodizeOrSmooth,ApodOrSmoothValue)
%
% Compute the noise of a signal (std of the signal). 
% - Do a linear fit of the sugnal and subtract that from signal. This is done because sometimes there is just 
%   simple linear baseline in the signal (eg for spectra if the water suppression did not work properly).
% - Divide the signal - LinearFit in several subregions, compute the stdev of all subregions and the stdev of the 
%   whole region.
% - Compute the ratio of the global and the local stddev. If this ratio is close to 1, there is just noise in the 
%   signal. Depending on this ratio, subtract a apodized or smoothed signal from the original signal. The 
%   smoothing/apodization is the stronger the closer the ratio is to 1, because if there is just noise inside
%   one only has to subtract something linear but no peaks.
%   With this subtraction the baseline of the signal is subtracted.
% - Compute the stddev of the Signal from which the smoothed/apodized signal was subtracted (so the baseline
%   corrected signal)
% Apodization: Filtering the fourier transformed signal (eg in the time domain) with an exponential function.

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    
ExtraSmoothPts = 10;
if(NoiseRegion(1) - 1 < ExtraSmoothPts || size(Signal,2) - NoiseRegion(2) < ExtraSmoothPts)             % if in the Signal there are not enough points to extend the NoiseRegion
   ExtraSmoothPts = min([NoiseRegion(1) - 1, size(Signal,2) - NoiseRegion(2)]);                         % then simply take as much as possible
end

PowerScalingValue_Smooth = 1.0;                                                                            % Value with which the ApodOrSmoothValue scales with the Noise_StdGlobalLocalRatio.
ConstScalingValue_Smooth = 8;
PowerScalingValue_Apod = 0.6;    %0.5                                                                         % Same for Apodization
ConstScalingValue_Apod = 16;      %4

% 0.3 Preparations

pause on

%% 1. Linear Fit & Subtraction


Polyfit_Signal_Real = polyfit(NoiseRegion(1)-ExtraSmoothPts:NoiseRegion(2)+ExtraSmoothPts,real(Signal(NoiseRegion(1)-ExtraSmoothPts:NoiseRegion(2)+ExtraSmoothPts)),1);
Polyfit_Signal_Imag = polyfit(NoiseRegion(1)-ExtraSmoothPts:NoiseRegion(2)+ExtraSmoothPts,imag(Signal(NoiseRegion(1)-ExtraSmoothPts:NoiseRegion(2)+ExtraSmoothPts)),1);

Signal_SubLinear = Signal - complex(Polyval(Polyfit_Signal_Real,1:size(Signal,2)),Polyval(Polyfit_Signal_Imag,1:size(Signal,2)));



%% 1. Compute Local Std
    
    % Compute the Ratio of the Std of the Whole Noise Region and Subregions of the Region
    TotalPeakregionPoints = NoiseRegion(2) - NoiseRegion(1) + 1;
    PointsPerSubregion = 40;
    TotalSubregions = floor(TotalPeakregionPoints / PointsPerSubregion);
    Noise_Subregions_Std = zeros([1 TotalSubregions]);
    
    

    for Subregion = 1:TotalSubregions-1
        Noise_Subregions_Std(Subregion) = std(real(Signal_SubLinear(NoiseRegion(1) + (Subregion-1)*PointsPerSubregion : NoiseRegion(1) + Subregion*PointsPerSubregion)));
    end
    Noise_Subregions_Std(Subregion + 1) = std(real(Signal_SubLinear(NoiseRegion(1) + Subregion*PointsPerSubregion : NoiseRegion(2))));    
    Noise_StdLocal = mean(Noise_Subregions_Std);
    
    
    
    %% 2. Exclude subregions with high std
    % Compute array that contains the indices of the excluded Subregions (eg [1 4]: subregions 1 and 4 gets excluded)
    ExcludeSubregions_Logical = Noise_Subregions_Std > Noise_StdLocal + 1.3*std(Noise_Subregions_Std);   % logical array indexing excluded subregion
    ExcludeSubregions_Index = find(ExcludeSubregions_Logical);                                         % Index of excluded subregion
    IncludeSubregions_Index = setdiff(1:TotalSubregions,ExcludeSubregions_Index);
        
    % Compute the Starting and end points of the excluded regions
    ExcludeSubregions_StartEndPoints = cell([1 sum(ExcludeSubregions_Logical)]);
    ExcludeSubregions_Points = [];
    for ExcludeSubregions_LoopIndex = 1:sum(ExcludeSubregions_Logical)
        if(ExcludeSubregions_Index(ExcludeSubregions_LoopIndex) == TotalSubregions)
            ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex} = [NoiseRegion(1) + (ExcludeSubregions_Index(ExcludeSubregions_LoopIndex)-1) * PointsPerSubregion,NoiseRegion(2)];
        else
            ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex} = [NoiseRegion(1) + (ExcludeSubregions_Index(ExcludeSubregions_LoopIndex)-1) * PointsPerSubregion, ...
                                                                             NoiseRegion(1) + ExcludeSubregions_Index(ExcludeSubregions_LoopIndex) * PointsPerSubregion - 1];
        end
    % Compute the vector indexing all excluded points, and that indexing all included
    ExcludeSubregions_Points = union(ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(1):ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(2),ExcludeSubregions_Points);
    end
    IncludeSubregions_Points = setdiff(NoiseRegion(1):NoiseRegion(2),ExcludeSubregions_Points);

    
%% 3. Compute the Apodization & Smoothing Values

if(~exist('ApodOrSmoothValue','var') || ApodOrSmoothValue == 0)    
    % Compute again the local stddev and the global std, both without the excluded subregions
    Noise_StdLocal = mean(Noise_Subregions_Std(IncludeSubregions_Index));
    Noise_StdGlobal = std(real(Signal_SubLinear(IncludeSubregions_Points)));
    Noise_StdGlobalLocalRatio = Noise_StdGlobal / Noise_StdLocal;



    if(strcmpi(ApodizeOrSmooth,'Apodize'))
        ApodOrSmoothValue = ConstScalingValue_Apod./abs(Noise_StdGlobalLocalRatio - 1).^PowerScalingValue_Apod;
    else
        ApodOrSmoothValue = round(ConstScalingValue_Smooth./abs(Noise_StdGlobalLocalRatio - 1).^PowerScalingValue_Smooth) + 10;
    end
    
end



%% 4. Smooth or apodize the data

if(strcmpi(ApodizeOrSmooth,'Apodize'))
    Signal_Time = ifft(ifftshift(Signal));
    Signal_Time_Apod = Signal_Time .* exp(-ApodOrSmoothValue*((1:numel(Signal_Time))-1)/numel(Signal_Time));
    Signal_ApodSmooth = fftshift(fft(Signal_Time_Apod));
else
    Signal_ApodSmooth_Real = transpose(smooth(real(Signal_SubLinear(NoiseRegion(1)-ExtraSmoothPts:NoiseRegion(2)+ExtraSmoothPts)),ApodOrSmoothValue));
    Signal_ApodSmooth_Imag = transpose(smooth(imag(Signal_SubLinear(NoiseRegion(1)-ExtraSmoothPts:NoiseRegion(2)+ExtraSmoothPts)),ApodOrSmoothValue));
    Signal_ApodSmooth = Signal;
    Signal_ApodSmooth(NoiseRegion(1):NoiseRegion(2)) = complex(Signal_ApodSmooth_Real(ExtraSmoothPts + 1 : end - ExtraSmoothPts), Signal_ApodSmooth_Imag(ExtraSmoothPts + 1 : end - ExtraSmoothPts)) + ...    % Smoothing
                                                       complex(Polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)), Polyval(Polyfit_Signal_Imag,NoiseRegion(1):NoiseRegion(2)));                       % FitLinear
end



%% 5. Subtract the Smoothed or Apodized Data

Signal_SubPeak = Signal - Signal_ApodSmooth;


  %% 2. Exclude subregions with high std
%     % Compute array that contains the indices of the excluded Subregions (eg [1 4]: subregions 1 and 4 gets excluded)
%     
%     
%     
%     ExcludeSubregions_Logical = Noise_Subregions_Std > Noise_StdLocal + 2*std(Noise_Subregions_Std);   % logical array indexing excluded subregion
%     ExcludeSubregions_Index = find(ExcludeSubregions_Logical);                                         % Index of excluded subregion
%     IncludeSubregions_Index = setdiff(1:TotalSubregions,ExcludeSubregions_Index);
%         
%     % Compute the Starting and end points of the excluded regions
%     ExcludeSubregions_StartEndPoints = cell([1 sum(ExcludeSubregions_Logical)]);
%     ExcludeSubregions_Points = [];
%     for ExcludeSubregions_LoopIndex = 1:sum(ExcludeSubregions_Logical)
%         if(ExcludeSubregions_Index(ExcludeSubregions_LoopIndex) == TotalSubregions)
%             ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex} = [NoiseRegion(1) + (ExcludeSubregions_Index(ExcludeSubregions_LoopIndex)-1) * PointsPerSubregion,NoiseRegion(2)];
%         else
%             ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex} = [NoiseRegion(1) + (ExcludeSubregions_Index(ExcludeSubregions_LoopIndex)-1) * PointsPerSubregion, ...
%                                                                              NoiseRegion(1) + ExcludeSubregions_Index(ExcludeSubregions_LoopIndex) * PointsPerSubregion - 1];
%         end
%     % Compute the vector indexing all excluded points, and that indexing all included
%     ExcludeSubregions_Points = union(ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(1):ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(2),ExcludeSubregions_Points);
%     end
%     IncludeSubregions_Points = setdiff(NoiseRegion(1):NoiseRegion(2),ExcludeSubregions_Points);
% 
% 




%% 6. Compute the noise of the SubPeak data with PeakRegions excluded.

Noise_std = std(horzcat([real(Signal_SubPeak(IncludeSubregions_Points)),imag(Signal_SubPeak(IncludeSubregions_Points))]));



%% DEBUG MODE: plot

% NoiseLinFit_fig = figure;
% plot(real(Signal(NoiseRegion(1):NoiseRegion(2))))
% hold on
% plot(Polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)),'r')
% % plot X-es in Regions that are excluded
% % for ExcludeSubregions_LoopIndex = 1:sum(ExcludeSubregions_Logical)
% %     StartPoint = ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(1);
% %     EndPoint = ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(2);
% %     plot([StartPoint,EndPoint],[min(Signal(StartPoint:EndPoint)),max(Signal(StartPoint:EndPoint))],'g','Linewidth',1.8) % plot line increasing diagonal from left to right
% %     plot([StartPoint,EndPoint],[max(Signal(StartPoint:EndPoint)),min(Signal(StartPoint:EndPoint))],'g','Linewidth',1.8) % plot line decreasing diagonal from left to right
% % end
% hold off
% 
% 
% NoiseSubLinear_fig = figure;
% plot(real(Signal_SubLinear(NoiseRegion(1):NoiseRegion(2))))
% hold on
% % plot X-es in Regions that are excluded
% for ExcludeSubregions_LoopIndex = 1:sum(ExcludeSubregions_Logical)
%     StartPoint = ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(1);
%     EndPoint = ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(2);    
%     plot([StartPoint - NoiseRegion(1),EndPoint - NoiseRegion(1)],[min(real(Signal_SubLinear(StartPoint:EndPoint))),max(real(Signal_SubLinear(StartPoint:EndPoint)))],'r','Linewidth',2.0) % plot line increasing diagonal from left to right
%     plot([StartPoint - NoiseRegion(1),EndPoint - NoiseRegion(1)],[max(real(Signal_SubLinear(StartPoint:EndPoint))),min(real(Signal_SubLinear(StartPoint:EndPoint)))],'r','Linewidth',2.0) % plot line decreasing diagonal from left to right
% end
% hold off
% 
% 
% 
% 
% NoiseApodSmooth_fig = figure;
% if(strcmpi(ApodizeOrSmooth,'Apodize'))
%     plot(real(Signal(NoiseRegion(1):NoiseRegion(2))))
% else
%     plot(real(Signal_SubLinear))
% end
% hold on
% plot(real(Signal_ApodSmooth(NoiseRegion(1):NoiseRegion(2))),'g')
% hold off
% 
% SubPeak_fig = figure;
% plot(real(Signal_SubPeak(NoiseRegion(1):NoiseRegion(2))))
% 
% pause
% waitforbuttonpress
% close(NoiseSubLinear_fig,NoiseLinFit_fig,NoiseApodSmooth_fig,SubPeak_fig);



%% 7. THE END

pause off







