function ShowSequenceArbGrads(PathToTxtFile,inputvar2)
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
GradValues.data = reshape(GradValues.data,[size_Data(1) numel(GradValues.data)/size_Data(1)]);





%% 2. Plot gradient values

% Time Plot

% Open Figure
TimePlotFig = figure;

% Set Axis
MaxValues = max(abs(GradValues.data));
%Plot_Lims = [min(min(GradValues.data(:,1)),-0.1*MaxValues(1)) max(max(GradValues.data(:,1)),0.1*MaxValues(1))];
%Plot_Lims = cat(2,Plot_Lims,[min(min(GradValues.data(:,2)),-0.1*MaxValues(2)) max(max(GradValues.data(:,2)),0.1*MaxValues(2))]);
Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
axis(Plot_Lims)

% Labels on axes
xlabel('GradStrength x [mT/m]')
ylabel('GradStrength y [mT/m]')



% Plot in a Loop
hold on
pause on
for i = 1:size(GradValues.data,1)
    scatter(GradValues.data(i,1),GradValues.data(i,2),'filled','b')
    pause(1)
end
hold off
pause off



%% 3. Postparations

close all;
% fclose(fid)






