function QualityMeasure = kSpace_DistBtwMeasPts(PtsMeas, CellSize)
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



%% 0. PREPARATIONS

% Assign standard values to variables if nothing is passed to function.
if(nargin < 2)
    display([ char(10) 'Gimme more input, Ma''am!' char(10) ])
    return;
end 






%% 2. Replicate the Measured in such a way to get all combinations for the distances between them
% (for every non-measured point we want the distance to all measured ones)

PtsMeas_rep = repmat(transpose(PtsMeas),[1 numel(PtsMeas)]);
PtsMeas2_rep = repmat(PtsMeas,[numel(PtsMeas) 1]);




%% 3. Find the x-y coordinates of the measured and not measured points within the cell.

[PtsMeas_x, PtsMeas_y] = ind2sub(CellSize,PtsMeas_rep);
[PtsMeas2_x, PtsMeas2_y] = ind2sub(CellSize,PtsMeas2_rep);





%% 4. Compute the x and y distances

Dist_x = abs(PtsMeas_x - PtsMeas2_x);
Dist_y = abs(PtsMeas_y - PtsMeas2_y);




%% 5. Correct the distance.
% (Sometimes a certain non-measured point is close to a measured point that is in the neighbouring replicated cell (the real undersampling pattern is a replicated form of the elementary cell)
% E.g.:
%
%   Elementary Cell:
%   1 o o x
%   o 2 o o
%
%   The right upper not-measured point (x) is closer to the 1 (denoted below as 3) of a replicated point:
%   1 o o o | 1 o o o | 1 o o o
%   o 2 o o | o 2 o o | o 2 o o
%   ---------------------------
%   1 o o o | 1 o o x | 3 o o o
%   o 2 o o | o 2 o o | o 2 o o
%   ---------------------------
%   1 o o o | 1 o o o | 1 o o o
%   o 2 o o | o 2 o o | o 2 o o




MaxDistance_x = floor(CellSize(1)/2);   % It is not possible that two points are farther apart from each other in the replicated cell
MaxDistance_y = floor(CellSize(2)/2);

Dist_x_mod = mod(CellSize(1),Dist_x);
Dist_y_mod = mod(CellSize(1),Dist_y);

Dist_x(Dist_x > MaxDistance_x) = Dist_x_mod(Dist_x > MaxDistance_x);
Dist_y(Dist_y > MaxDistance_y) = Dist_y_mod(Dist_y > MaxDistance_y);





%% 6. Compute the Euclidean distance and the QualityMeasures

Distance_mat = sqrt(Dist_x.^2 + Dist_y.^2);
Distance_mat = triu(Distance_mat,1);        % Only take upper right part (without diagonal, therefore the 1) of the matrix
Distance_mat(Distance_mat == 0) = NaN;      % Distances of 0 should not count.




%% 7. Compute all QualityMeasures


% The minimum distance between points. 
QualityMeasure_min = nanmin(reshape(Distance_mat, [1 numel(Distance_mat)]));

% The mean distance between points.
QualityMeasure_std = nanstd(reshape(Distance_mat, [1 numel(Distance_mat)]));


QualityMeasure = [QualityMeasure_std QualityMeasure_min QualityMeasure_std/(QualityMeasure_min^2)];
%QualityMeasure = QualityMeasure_std/(QualityMeasure_min^1.5);

