function [GradValues, GradPos, GradSlew] = ShowSequenceArbGrads_mediumnew(PathToTxtFile)
%
% ShowSequenceArbGrads Do nothing specific
%
% This function was written by Bernhard Strasser, Jannuary 2015.
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         PathToTxtFile                   ...    Path to text file containing the arbitrary gradient values
%
% Output:
% -         A                           ...     This is output A
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('PathToTxtFile','var'))
    fprintf('\nNo input file. Abort.')
    return
end



% 0.3 Definitions
    





%% 1. Read File

GradValues = importdata(PathToTxtFile,'\t');


size_Data = size(GradValues.data);
GradValues.data(isnan(GradValues.data)) = [];
if(mod(numel(GradValues.data),3) == 1)
    GradValues.data(end) = [];
    size_Data(1) = size_Data(1)-1;
end
GradValues.data = reshape(GradValues.data,[size_Data(1) size_Data(2)-2]);


%% 2. Calculate Positions

% Need to integrate, i.e. s(t_n) = sum(i<n, v(t_i) * Delta_t) = Delta_t * sum(i<n,v(t_i))            % really i<n or i<=n ?
GradPos_x = zeros([1 size_Data(1)]);
GradPos_y = zeros([1 size_Data(1)]);

try
    GradPos_x = sum(triu(repmat(GradValues.data(:,1),[size_Data(1) size_Data(1)])));
    GradPos_y = sum(triu(repmat(GradValues.data(:,2),[size_Data(1) size_Data(1)])));
catch errie
    GradPos_x(1) = GradValues.data(1,1);
    GradPos_y(1) = GradValues.data(1,2);
    for i=2:size_Data(1)
        GradPos_x(i) = GradPos_x(i-1) + GradValues.data(i,1);
        GradPos_y(i) = GradPos_y(i-1) + GradValues.data(i,2);        
    end
end


GradPos = transpose(cat(2,[0; 0],cat(1,GradPos_x,GradPos_y))) * 10;     % in units of mT*us/m
clear GradPos_x GradPos_y




%% 3. Calculate SlewRates

% Numeric differentiation: a(t_n) = (a(t_n) - a(t_n-1)) / Delta_t

GradSlew = GradValues.data - circshift(GradValues.data,[-1 0]);        % Then the last value is not meaningful
GradSlew(end,:) = [];
GradSlew = GradSlew(:,1:2) / 10;                                              % in units of mT/m per us



%% 3. Plot gradient values

% Time Plot

% Grad Values
GradFig = figure;
hold on
movegui(GradFig,'northwest')

% Set Axis
MaxValues = max(abs(GradValues.data));
Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
axis(Plot_Lims)

% Labels on axes
xlabel('GradStrength x [mT/m]')
ylabel('GradStrength y [mT/m]')



% Grad Pos
GradPosFig = figure;
hold on
movegui(GradPosFig,'northeast')

% Set Axis
MaxValues = max(abs(GradPos));
Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
axis(Plot_Lims)

% Labels on axes
xlabel('GradPos x [mT*us/m]')
ylabel('GradPos y [mT*us/m]')




% Grad SlewRate
GradSlewFig = figure;
hold on
movegui(GradSlewFig,'south')

% Set Axis
MaxValues = max(abs(GradSlew));
Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
axis(Plot_Lims)

% Labels on axes
xlabel('GradSlew x [mT/m per us]')
ylabel('GradSlew y [mT/m per us]')




% Plot in a Loop
pause on
for i = 1:size(GradPos,1)
    if(i>1)
        figure(GradFig)
        scatter(GradValues.data(i-1,1),GradValues.data(i-1,2),'filled','b')
    end

    figure(GradPosFig)
    scatter(GradPos(i,1),GradPos(i,2),'filled','b')
    
    if(i>2)
        figure(GradSlewFig)
        scatter(GradSlew(i-2,1),GradSlew(i-2,2),'filled','b')  
        sqrt(GradSlew(i-2,:) * transpose(GradSlew(i-2,:)))
    end
    
    pause
end
hold off
pause off



%% 3. Postparations

% fclose(fid)






