function [BestPatterns,no_Patterns,QualityMeasure,R_InPlaneVD] = FindBestGrappaPattern(R_InPlane)
%
% kSpace_Distance Compute quality measure of kSpace Pattern.
%
% This function was written by Bernhard Strasser, July 2013.
%
%
% This function 1) calculates the distance of a non-measured kSpace point to all the measured ones.
%               2) Does this for all non-measured points.
%               3) Calculates a quality measure for those distances (like maximum, mean or combination).
%               4) Does this for all non-measured points.
%               5) Computes the same quality measure using all non-measured points (maximum, mean or combination).
%
%
% [csi,csi_kspace] = read_csi_1_4(csi_path, zerofill_to_nextpow2_flag, zerofilling_fact, x_shift,y_shift)
%
% Input: 
% -         PtsMeas                     ...     The measured points as linear index (e.g. if A=[1 0;0 1], 1...measured --> PtsMeas = [1 4])
% -         CellSize                    ...     The size of the elementary caipirinha cell. in the above example: CellSize = [2 2]
%                                               Use at least 2*[Rx,Ry], where Rx and Ry would be the GRAPPA acceleration factors in x and y dir.
%                                               Otherwise you get only 1 possible undersampling cell.
%
% Output:
% -         QualityMeasure              ...     The quality measure of the kSpace Pattern, e.g. the mean distance of non-measured to measured points,
%                                               or the maximum of the distances between measured and non-measured ones.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: none yet.




%% 0. DEFINITIONS, PREPARATIONS





%% 1. Create All Possible Patterns of the Measured Points.

[R_InPlane_x, R_InPlane_y] = meshgrid(1:4);

R_InPlane_x = reshape(R_InPlane_x,[1 numel(R_InPlane_x)]);
R_InPlane_y = reshape(R_InPlane_y,[1 numel(R_InPlane_y)]);




%% 2. Restrict Patterns

% Criteria: - If one acceleration is more than four times the other --> dont allow
%			- If the R_InPlane is not reached --> dont allow
%			- If R_InPlaneActual is much higher than R_InPlane





% Omit too big differences between R_x and R_x
XDoubleY = R_InPlane_x > 4*R_InPlane_y;
YDoubleX = R_InPlane_y > 4*R_InPlane_x;

R_InPlane_Actual = R_InPlane_x .* R_InPlane_y;

% If the R_InPlane is not reached --> dont allow
R_InPlaneActualTooSmall = (R_InPlane_Actual < R_InPlane);

% If R_InPlaneActual is much higher than R_InPlane
R_InPlaneActualTooLarge = (R_InPlane_Actual > 1.5*R_InPlane);




R_InPlane_x(XDoubleY | YDoubleX | R_InPlaneActualTooSmall | R_InPlaneActualTooLarge) = [];
R_InPlane_y(XDoubleY | YDoubleX | R_InPlaneActualTooSmall | R_InPlaneActualTooLarge) = [];
R_InPlane_Actual(XDoubleY | YDoubleX | R_InPlaneActualTooSmall | R_InPlaneActualTooLarge) = [];

R_InPlaneVD = 1./(1/R_InPlane - 1./R_InPlane_Actual); 


%% 3. Define Important Variables and Quality Measure

BestPatterns.Indices = 1:numel(R_InPlane_x);
BestPatterns.Patterns = ones([1 numel(R_InPlane_x)]);
BestPatterns.CellSizes = cat(2,R_InPlane_x,R_InPlane_y);
no_Patterns = numel(R_InPlane_x);



QltyMeas_dummy = DistBtwMeasPts(squeeze(BestPatterns.Patterns(1,:)), BestPatterns.CellSizes,0);
NoQltyMeas = size(QltyMeas_dummy,2);
QualityMeasure = zeros([size(BestPatterns.Patterns,1) NoQltyMeas]);
clear QltyMeas_dummy NoQltyMeas

for Patt_no = 1:size(BestPatterns.Patterns,1)
    
    QualityMeasure(Patt_no,:) = DistBtwMeasPts(squeeze(BestPatterns.Patterns(Patt_no,:)), BestPatterns.CellSizes,0);
    
end









