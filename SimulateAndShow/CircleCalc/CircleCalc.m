function CircleCoordinates = CircleCalc(Radius,CircleCenter,NoOfPointsToCalculate)
%
% CircleCalc Calculate coordinated of a circle
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% .
% 
%
%
% CircleCoordinates = CircleCalc(Radius,CircleCenter,NoOfPointsToCalculate)
%
% Input: 
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         A                           ...     This is output A
% -         B                           ...     This is output B
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('Radius','var'))
    Radius = 1;
end
if(~exist('CircleCenter','var'))
    CircleCenter = [0,0];
end
if(~exist('NoOfPointsToCalculate','var'))
    NoOfPointsToCalculate = 100;
end
if(NoOfPointsToCalculate > 10^5)
    fprintf('\nWarning: Too many points for drawing circle demanded. Will restrict to 10^5.')
    NoOfPointsToCalculate = 10^5;
end

% 0.3 Definitions
    






%% 1. Calculate Angle-Vector 

% We need exactly one turn, so lets go from 0 till 2*pi-Fineity, where Fineity is the angle between two points on the circle
Fineity = 2*pi/(NoOfPointsToCalculate-1);
phi = 0:Fineity:(2*pi);


%% 2. Calculate circle

CircleCoordinates = Radius*[sin(phi);cos(phi)] + repmat(reshape(CircleCenter,[2 1]),[1 NoOfPointsToCalculate]);




%% 3. Postparations

% fclose(fid)






