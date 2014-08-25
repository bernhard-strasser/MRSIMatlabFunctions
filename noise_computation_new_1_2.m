function [Noise_std,Signal_ApodSmooth,Signal_SubSmooth,ExcludeSubregions_StartEndPoints] = noise_computation_new_1_2(Signal,NoiseRegion,ApodizeOrSmooth,ApodOrSmoothValue)
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
    
PowerScalingValue_Smooth = 0.5;  %0.5                                                                            % Value with which the ApodOrSmoothValue scales with the Noise_StdGlobalLocalRatio.
ConstScalingValue_Smooth = 3.0;  %3.0
PowerScalingValue_Apod = 0.6;    %0.5   0.6      0.47   0.6                                                             % Same for Apodization
ConstScalingValue_Apod = 6.9;      %4    16       3.2   6.9

% 0.3 Preparations

NoOfNoiseRegions = size(NoiseRegion,2);
IncludeSubregions_AllPoints = [];

ExcludeSubregions_StartEndPoints = cell(1,NoOfNoiseRegions);
Signal_ApodSmooth = cell(1,NoOfNoiseRegions);
Signal_SubSmooth = cell(1,NoOfNoiseRegions);


ExtraSmoothPts = cell(1,NoOfNoiseRegions);
for NoiseIndex = 1:NoOfNoiseRegions
    ExtraSmoothPts{NoiseIndex} = 10;
    if(NoiseRegion{NoiseIndex}(1) - 1 < ExtraSmoothPts{NoiseIndex} || size(Signal,2) - NoiseRegion{NoiseIndex}(2) < ExtraSmoothPts{NoiseIndex})     % if in the Signal there are not enough points to extend the NoiseRegion
       ExtraSmoothPts{NoiseIndex} = min([NoiseRegion{NoiseIndex}(1) - 1, size(Signal,2) - NoiseRegion{NoiseIndex}(2)]);                             % then simply take as much as possible
    else
    end
end

pause on



%% Iteration over all NoiseRegions

for NoiseIndex = 1:NoOfNoiseRegions

    %% 1. Linear Fit & Subtraction

    ExtraPts = ExtraSmoothPts{NoiseIndex};
    NoiseReg = NoiseRegion{NoiseIndex};

    Polyfit_Signal_Real = polyfit(NoiseReg(1)-ExtraPts:NoiseReg(2)+ExtraPts,real(Signal(NoiseReg(1)-ExtraPts:NoiseReg(2)+ExtraPts)),1);
    Polyfit_Signal_Imag = polyfit(NoiseReg(1)-ExtraPts:NoiseReg(2)+ExtraPts,imag(Signal(NoiseReg(1)-ExtraPts:NoiseReg(2)+ExtraPts)),1);

    Signal_SubLinear = Signal - complex(polyval(Polyfit_Signal_Real,1:size(Signal,2)),polyval(Polyfit_Signal_Imag,1:size(Signal,2)));


    %% 2. Compute Local Std

    % Compute the Ratio of the Std of the Whole Noise Region and Subregions of the Region
    TotalPeakregionPoints = NoiseReg(2) - NoiseReg(1) + 1;
    PointsPerSubregion = 25;
    TotalSubregions = floor(TotalPeakregionPoints / PointsPerSubregion);
    Noise_Subregions_Std = zeros([1 TotalSubregions]);
    %Noise_Subregions_Mean = zeros([1 TotalSubregions]);    


    for Subregion = 1:TotalSubregions-1
        Noise_Subregions_Std(Subregion) = std(real(Signal_SubLinear(NoiseReg(1) + (Subregion-1)*PointsPerSubregion : NoiseReg(1) + Subregion*PointsPerSubregion)));
        %Noise_Subregions_Mean(Subregion) = mean(real(Signal_SubLinear(NoiseReg(1) + (Subregion-1)*PointsPerSubregion : NoiseReg(1) + Subregion*PointsPerSubregion)));
    end
    Noise_Subregions_Std(Subregion + 1) = std(real(Signal_SubLinear(NoiseReg(1) + Subregion*PointsPerSubregion : NoiseReg(2))));
    %Noise_Subregions_Mean(Subregion + 1) = mean(real(Signal_SubLinear(NoiseReg(1) + Subregion*PointsPerSubregion : NoiseReg(2))));   
    Noise_StdLocal = mean(Noise_Subregions_Std);





    %% 3. Compute the Apodization & Smoothing Values

    if(~exist('ApodOrSmoothValue','var') || ApodOrSmoothValue == 0)    
        % Compute again the local stddev and the global std, both without the excluded subregions
        Noise_StdGlobal = std(real(Signal_SubLinear(NoiseReg(1):NoiseReg(2))));
        Noise_StdGlobalLocalRatio = Noise_StdGlobal / Noise_StdLocal;


        if(strcmpi(ApodizeOrSmooth,'Apodize'))
            ApodOrSmoothValue = ConstScalingValue_Apod./abs(Noise_StdGlobalLocalRatio - 1).^PowerScalingValue_Apod + 25
        else
            ApodOrSmoothValue = round(ConstScalingValue_Smooth./abs(Noise_StdGlobalLocalRatio - 1).^PowerScalingValue_Smooth) + 58;
        end

    end



    %% 4. Smooth or apodize the data

    if(strcmpi(ApodizeOrSmooth,'Apodize'))
        Signal_Time = ifft(ifftshift(Signal));
        Signal_Time_Apod = Signal_Time .* exp(-ApodOrSmoothValue*((1:numel(Signal_Time))-1)/numel(Signal_Time));
        Signal_ApodSmooth{NoiseIndex} = fftshift(fft(Signal_Time_Apod));
    else
        Signal_ApodSmooth_Real = transpose(smooth(real(Signal_SubLinear(NoiseReg(1)-ExtraPts:NoiseReg(2)+ExtraPts)),ApodOrSmoothValue));
        Signal_ApodSmooth_Imag = transpose(smooth(imag(Signal_SubLinear(NoiseReg(1)-ExtraPts:NoiseReg(2)+ExtraPts)),ApodOrSmoothValue));
        Signal_ApodSmooth{NoiseIndex} = Signal;
        Signal_ApodSmooth{NoiseIndex}(NoiseReg(1):NoiseReg(2)) = ...
            complex(Signal_ApodSmooth_Real(ExtraPts + 1 : end - ExtraPts), Signal_ApodSmooth_Imag(ExtraPts + 1 : end - ExtraPts)) + ...                                     % Smoothing
            complex(polyval(Polyfit_Signal_Real,NoiseReg(1):NoiseReg(2)), polyval(Polyfit_Signal_Imag,NoiseReg(1):NoiseReg(2)));                % FitLinear
    end

    Signal_SubSmooth{NoiseIndex} = Signal - Signal_ApodSmooth{NoiseIndex};




    %% 5. Exclude subregions with high std
    % Compute array that contains the indices of the excluded Subregions (eg [1 4]: subregions 1 and 4 gets excluded)

    % for Subregion = 1:TotalSubregions-1
    %     Noise_Subregions_Std(Subregion) = std(real(Signal_SubSmooth(NoiseRegion(1) + (Subregion-1)*PointsPerSubregion : NoiseRegion(1) + Subregion*PointsPerSubregion)));
    %     Noise_Subregions_Mean(Subregion) = mean(real(Signal_SubSmooth(NoiseRegion(1) + (Subregion-1)*PointsPerSubregion : NoiseRegion(1) + Subregion*PointsPerSubregion)));
    % end
    % Noise_Subregions_Std(Subregion + 1) = std(real(Signal_SubSmooth(NoiseRegion(1) + Subregion*PointsPerSubregion : NoiseRegion(2))));
    % Noise_Subregions_Mean(Subregion + 1) = mean(real(Signal_SubSmooth(NoiseRegion(1) + Subregion*PointsPerSubregion : NoiseRegion(2))));   
    % Noise_StdLocal = mean(Noise_Subregions_Std);
    % 
    % %ExcludeSubregions_Logical = abs(Noise_Subregions_Mean) > 0.25 * span(Signal_SubSmooth(NoiseRegion(1):NoiseRegion(2)));  
    % %ExcludeSubregions_Logical = Noise_Subregions_Mean > mean(Noise_Subregions_Mean) + 1.4*std(Noise_Subregions_Mean);   
    % %ExcludeSubregions_Logical = Noise_Subregions_Std > Noise_StdLocal + 1.6*std(Noise_Subregions_Std);
    % ExcludeSubregions_Logical = Noise_Subregions_Std > Noise_StdLocal + 1.8*std(Noise_Subregions_Std); %| ...   % logical array indexing excluded subregion; Exclude regions that deviate more than 1.5*std from the mean.
    %                             %Noise_Subregions_Mean > mean(Noise_Subregions_Mean) + 1.4*std(Noise_Subregions_Mean);  
    % ExcludeSubregions_Index = find(ExcludeSubregions_Logical);                                                % Index of excluded subregion. within (mean +- 1*std) 68.27% of all values are inside, 
    % IncludeSubregions_Index = setdiff(1:TotalSubregions,ExcludeSubregions_Index);                             % within +- 2*std there are 95.45% of all values, if it is gaussian distributed.
    % 
    % % Compute the Starting and end points of the excluded regions
    % ExcludeSubregions_StartEndPoints = cell([1 sum(ExcludeSubregions_Logical)]);
    % ExcludeSubregions_Points = [];
    % for ExcludeSubregions_LoopIndex = 1:sum(ExcludeSubregions_Logical)
    %     if(ExcludeSubregions_Index(ExcludeSubregions_LoopIndex) == TotalSubregions)
    %         ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex} = [NoiseRegion(1) + (ExcludeSubregions_Index(ExcludeSubregions_LoopIndex)-1) * PointsPerSubregion,NoiseRegion(2)];
    %     else
    %         ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex} = [NoiseRegion(1) + (ExcludeSubregions_Index(ExcludeSubregions_LoopIndex)-1) * PointsPerSubregion, ...
    %                                                                          NoiseRegion(1) + ExcludeSubregions_Index(ExcludeSubregions_LoopIndex) * PointsPerSubregion - 1];
    %     end
    %   % Compute the vector indexing all excluded points, and that indexing all included
    %   ExcludeSubregions_Points = union(ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(1):ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(2),ExcludeSubregions_Points);
    % end
    % IncludeSubregions_Points = setdiff(NoiseRegion(1):NoiseRegion(2),ExcludeSubregions_Points);


    bla_smooth = 0;
    bla_smooth(NoiseReg(1):NoiseReg(2)) = transpose(SmoothAlsoEndPts_1_0(real(Signal_SubSmooth{NoiseIndex}(NoiseReg(1):NoiseReg(2))),14));
    PotentialExcludePts_Logical = abs(bla_smooth) > 0.6 * std(Signal_SubSmooth{NoiseIndex}(NoiseReg(1):NoiseReg(2)));
    ExcludeSubregions_StartEndPoints{NoiseIndex} = findgroup_1_1(PotentialExcludePts_Logical,10);
    ExcludeSubregions_Points = [];
    for Loopy = 1:size(ExcludeSubregions_StartEndPoints{NoiseIndex},2)
        ExcludeSubregions_Points = union(ExcludeSubregions_Points, ExcludeSubregions_StartEndPoints{NoiseIndex}{Loopy}(1) : ExcludeSubregions_StartEndPoints{NoiseIndex}{Loopy}(2));
    end
    

    IncludeSubregions_Points = setdiff(NoiseReg(1):NoiseReg(2),ExcludeSubregions_Points);
    IncludeSubregions_AllPoints = union(IncludeSubregions_AllPoints,IncludeSubregions_Points);


    % std(Signal_SubSmooth(NoiseRegion(1):NoiseRegion(2))
    % bla_smooth = transpose(smooth(real(Signal_SubSmooth(NoiseRegion(1):NoiseRegion(2))),14));
    % PotentialExcludePts_Logical = abs(bla_smooth) > 0.5 * std(Signal_SubSmooth(NoiseRegion(1):NoiseRegion(2)));
    % ExcludeSubregions_StartEndPoints = findgroup_1_1(PotentialExcludePts_Logical,12);
    % ExcludeSubregions_Points = [];
    % for Loopy = 1:size(ExcludeSubregions_StartEndPoints,2)
    %     ExcludeSubregions_StartEndPoints{Loopy} = ExcludeSubregions_StartEndPoints{Loopy} + NoiseRegion(1) - 1;
    %     ExcludeSubregions_Points = union(ExcludeSubregions_Points, ExcludeSubregions_StartEndPoints{Loopy}(1) : ExcludeSubregions_StartEndPoints{Loopy}(2));
    % end
    % 
    % IncludeSubregions_Points = setdiff(NoiseRegion(1):NoiseRegion(2),ExcludeSubregions_Points);


end



%% 6. Compute the noise of the SubPeak data with PeakRegions excluded.

% Concatenate all the Noises


Noise_std = std(horzcat([real(Signal_SubSmooth{NoiseIndex}(IncludeSubregions_AllPoints)),imag(Signal_SubSmooth{NoiseIndex}(IncludeSubregions_AllPoints))]));


%% DEBUG MODE: plot

% NoiseLinFit_fig = figure;
% plot(real(Signal(NoiseRegion(1):NoiseRegion(2))))
% hold on
% plot(polyval(Polyfit_Signal_Real,NoiseRegion(1):NoiseRegion(2)),'r')
% hold off
% 
% 
% NoiseSubLinear_fig = figure;
% plot(real(Signal_SubLinear(NoiseRegion(1):NoiseRegion(2))))
% 
% 
% NoiseApodSmooth_fig = figure;
% if(strcmpi(ApodizeOrSmooth,'Apodize'))
%     plot(real(Signal(NoiseRegion(1):NoiseRegion(2))))
% else
%     plot(real(Signal_SubLinear(NoiseRegion(1):NoiseRegion(2))))
% end
% hold on
% plot(real(Signal_ApodSmooth(NoiseRegion(1):NoiseRegion(2))),'r')
% hold off
% 
% SubSmooth_fig = figure;
% plot(real(Signal_SubSmooth(NoiseRegion(1):NoiseRegion(2))))
% 
% 
% 
% 
% SubPeak_fig = figure;
% plot(real(Signal_SubSmooth(NoiseRegion(1):NoiseRegion(2))))
% hold on
% plot(abs(bla_smooth(NoiseRegion(1):NoiseRegion(2))),'g')
% % plot X-es in Regions that are excluded
% for ExcludeSubregions_LoopIndex = 1:size(ExcludeSubregions_StartEndPoints,2)
%     StartPoint = ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(1);
%     EndPoint = ExcludeSubregions_StartEndPoints{ExcludeSubregions_LoopIndex}(2);    
%     plot([StartPoint - NoiseRegion(1),EndPoint - NoiseRegion(1)],[min(real(Signal_SubLinear(StartPoint:EndPoint))),max(real(Signal_SubLinear(StartPoint:EndPoint)))],'r','Linewidth',2.0) % plot line increasing diagonal from left to right
%     plot([StartPoint - NoiseRegion(1),EndPoint - NoiseRegion(1)],[max(real(Signal_SubLinear(StartPoint:EndPoint))),min(real(Signal_SubLinear(StartPoint:EndPoint)))],'r','Linewidth',2.0) % plot line decreasing diagonal from left to right
% end
% hold off
% 
% 
% Noise_fig = figure;
% plot(horzcat([real(Signal_SubSmooth(IncludeSubregions_Points)),imag(Signal_SubSmooth(IncludeSubregions_Points))]))
% 
% 
% pause
% waitforbuttonpress
% close(NoiseSubLinear_fig,NoiseLinFit_fig,NoiseApodSmooth_fig,SubSmooth_fig,SubPeak_fig,Noise_fig);



%% 7. THE END

pause off







