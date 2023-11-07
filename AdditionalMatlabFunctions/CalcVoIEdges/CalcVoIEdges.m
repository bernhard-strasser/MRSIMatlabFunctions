function [VOI_X_MIN,VOI_X_MAX] = CalcVoIEdges(MatrixSize,VoISize,FoVSize,ExtraVoxels,RoundUpOrDown)
%
% CalcVoIEdges Calculate the edge voxels of the VoI
%
% This function was written by Bernhard Strasser, May 2018.
%
%
% The function can calculate the edge voxels of the VoI given the Matrix size, FoVSize and the number of extra voxels that
% should be included on each side of the VoI.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         MatrixSize                   ...    Self explaining, right?
% -         VoISize                      ...    Self explaining, right?
% -         FoVSize                      ...    Self explaining, right?
% -         ExtraVoxels                  ...    How many ExtraVoxels on each side of the VoI should be used

%
% Output:
% -         VOI_X_MIN                           ...     The calculated minimum voxel at the edge of the VoI (- ExtraVoxels)
% -         VOI_X_MAX                           ...     The calculated maximum voxel at the edge of the VoI (+ ExtraVoxels)
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

% if(~exist('inputvar1','var'))
%     inputvar1 = 0;
% end



% 0.3 Definitions
    

Oddiness = mod(MatrixSize,2);

if(~exist('RoundUpOrDown','var'))
    RoundUpOrDown = 'RoundUp';
end
if(strcmpi(RoundUpOrDown,'RoundDown'))
    RoundDownFlag = true;
else
    RoundDownFlag = false;    
end



%% 1. 


% If we have an even number of phase encoding steps, say 10, then the center of the VoI is in between two steps, in this case
% between step 5 and 6. If we have an VoI of only 1 step, we actually excited only half of slice 5 and half of 6, so nothing
% (but even in that case we want to restrict it to have 2 slices, otherwise we get errors)
% If we have an VoI of 2 slices, we fully encode the 2 slices. If we have an VoI of 3 slices, we encode just 2 slices plus 1/2
% slice on each side, so we just take the 2 full slices. 
% etc.
% For an odd number of phase encoding steps, the VoI-center is also at the FoV-center, and we can only have odd number of fully
% excited VoI-voxels.

% Let R be the number of times the voxel-size fits into the VoISize
% The following formula maps R to the NoOfVoxels as follows (for RoundDownFlag = false):
% Odd MatrixSizes:  0<=R<2.0 --> 0; 2.0<=R<4.0 --> 2; 4.0<=R<6.0 --> 4; ...
% Even MatrixSizes: 0<=R<1.0 --> 0; 1.0<=R<3.0 --> 2; 3.0<=R<5.0 --> 4; ...
NoOfVoxels = 2*floor(VoISize/(2*FoVSize/MatrixSize) + 0.5*(1-Oddiness) -double(RoundDownFlag)*0.5);

% Make NoOfVoxels odd for odd MatrixSizes
NoOfVoxels = NoOfVoxels + Oddiness;

% Add ExtraVoxels
NoOfVoxels = NoOfVoxels + 2*ExtraVoxels;

% Restrict NoOfVoxels to 1 for odd MatrixSizes, and 2 for even
NoOfVoxels(NoOfVoxels < (2-Oddiness)) = (2-Oddiness); 

% Restrict NoOfVoxels to MatrixSize
NoOfVoxels(NoOfVoxels > MatrixSize) = MatrixSize; 

% In total, we mapped R to the NoOfVoxels as follows for ExtraVoxels = 0:
% Odd MatrixSizes:  0<=R<2.0 --> 1; 2.0<=R<4.0 --> 3; 4.0<=R<6.0 --> 5; ...
% Even MatrixSizes: 0<=R<3.0 --> 2; 3.0<=R<5.0 --> 4; 5.0<=R<7.0 --> 6; ...


% Calculate the min and max voxel
Ctr = floor(MatrixSize/2)+1;
VOI_X_MIN = Ctr - floor(NoOfVoxels/2);
VOI_X_MAX = Ctr + ceil(NoOfVoxels/2) - 1;






%% 3. Postparations







