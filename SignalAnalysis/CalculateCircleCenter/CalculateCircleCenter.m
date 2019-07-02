function [CircleCenter,Radius] = CalculateCircleCenter(CirclePoints,BeginAndEndOfCircle)
%
% template_1_0 Do nothing specific
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
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

if(~exist('CirclePoints','var'))
    CircleCenter = 0; Radius = 0;
    fprintf('\nToo little input.')
    return
end

if(exist('BeginAndEndOfCircle','var'))
    if(BeginAndEndOfCircle(1) > 1)
        CirclePoints = CirclePoints(BeginAndEndOfCircle(1):end,:);
    end
    if(numel(BeginAndEndOfCircle) > 1 && BeginAndEndOfCircle(2) > 1)
        CirclePoints = CirclePoints(1:BeginAndEndOfCircle(2),:);
    end    
end
if(size(CirclePoints,1) < 3)
    CircleCenter = 0; Radius = 0;
    fprintf('\nToo few points to calculate circle center. Minimum 3.')
    return
end
if(size(CirclePoints,2) ~= 2)
    CircleCenter = 0; Radius = 0;
    fprintf('\nNeed x- and y- component of points in second dimension of CirclePoints.')
    return
end

% 0.3 Definitions
    






%% 1. Compute Center

% Take the following points: 1, 1/4*size, 3/5*size
Point1 = 1; 
Point2 = round(size(CirclePoints,1)/4)+1; 
Point3 = round(3*size(CirclePoints,1)/5)+1;


mr = (CirclePoints(Point2,2)-CirclePoints(Point1,2))/(CirclePoints(Point2,1)-CirclePoints(Point1,1));
mt = (CirclePoints(Point3,2)-CirclePoints(Point2,2))/(CirclePoints(Point3,1)-CirclePoints(Point2,1));
CircleCenter(1) = (mr*mt*(CirclePoints(Point3,2)-CirclePoints(Point1,2)) + mr*(CirclePoints(Point2,1)+CirclePoints(Point3,1)) - mt*(CirclePoints(Point1,1)+CirclePoints(Point2,1)))/(2*(mr-mt));
CircleCenter(2) = (CirclePoints(Point1,2)+CirclePoints(Point2,2))/2 - 1/mr*(CircleCenter(1)-(CirclePoints(Point1,1)+CirclePoints(Point2,1))/2);




%% 2. Compute Radius

Radius =  sqrt(sum((CirclePoints(end,:) - CircleCenter).^2));




%% 3. Postparations

% fclose(fid)






