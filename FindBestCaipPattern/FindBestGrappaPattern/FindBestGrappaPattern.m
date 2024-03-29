function [BestPatterns,no_Patterns,R_VD] = FindBestGrappaPattern(R,LookUpTable)
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



%% 1. LookUpTable

if(~exist('LookUpTable','var') || LookUpTable)
	
	if(ceil(R) == 2)
		R_x = [1 2];
		R_y = [2 1];
	elseif(ceil(R) == 3)
		R_x = [1 2 3];
		R_y = [3 2 1];	
	elseif(ceil(R) == 4)
		R_x = [1 2 4];
		R_y = [4 2 1];		
	elseif(ceil(R) == 5)
		R_x = [2 3];
		R_y = [3 2];
	elseif(ceil(R) == 6)
		R_x = [2 3];
		R_y = [3 2];	
	elseif(ceil(R) == 7)
		R_x = [2 3 4];
		R_y = [4 3 2];	
	elseif(ceil(R) == 8)
		R_x = [2 3 4];
		R_y = [4 3 2];	
	elseif(ceil(R) == 9)
		R_x = 3;
		R_y = 3;
	elseif(ceil(R) == 10)
		R_x = [3 4];
		R_y = [4 3];
	end
	
	R_Actual = R_x .* R_y;

	
end
	
	

%% 2. Routine


if(exist('LookUpTable','var') && ~LookUpTable)
	% 1. Create All Possible Patterns of the Measured Points.

	[R_x, R_y] = meshgrid(1:4);

	R_x = reshape(R_x,[1 numel(R_x)]);
	R_y = reshape(R_y,[1 numel(R_y)]);




	% 2. Restrict Patterns

	% Criteria: - If one acceleration is more than four times the other --> dont allow
	%			- If the R is not reached --> dont allow
	%			- If R_Actual is much higher than R_


	% Omit too big differences between R_x and R_x
	XDoubleY = R_x > 4*R_y;
	YDoubleX = R_y > 4*R_x;

	R_Actual = R_x .* R_y;

	% If the R is not reached --> dont allow
	R_ActualTooSmall = (R_Actual < R);

	% If R_Actual is much higher than R_
	R_ActualTooLarge = (R_Actual > 1.5*R);

	R_x(XDoubleY | YDoubleX | R_ActualTooSmall | R_ActualTooLarge) = [];
	R_y(XDoubleY | YDoubleX | R_ActualTooSmall | R_ActualTooLarge) = [];
	R_Actual(XDoubleY | YDoubleX | R_ActualTooSmall | R_ActualTooLarge) = [];

	
end


%% 3. Define Important Variables and Quality Measure


R_VD = 1./(1/R - 1./R_Actual); 


BestPatterns.Indices = 1:numel(R_x);
BestPatterns.Patterns = ones([1 numel(R_x)]);
BestPatterns.CellSizes = cat(1,R_x,R_y); BestPatterns.CellSizes = transpose(BestPatterns.CellSizes);
no_Patterns = numel(R_x);












