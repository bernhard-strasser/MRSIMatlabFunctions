function [OutData,FoV_Phase] = kSpace_FoVShift(InData, FoV_shifts) 
% 
% kSpace_FoVShift Perform In-Plane Shifts of the FoV of MR(S)I Data in kSpace
% 
%  [OutData,FoV_Phase] = kSpace_FoVShift(InData, FoV_shifts) 
%       
%   Input:      
% -     InData             Unshifted Data            (size: [#coils, nx, ny, nSlice, nTime]) (nTime = 1 for MRI)
%                          Must have as many slices as the OutData should have. Use similar sequence parameters as for the InData.
% -     FoV_shifts         The multiples of the FoV with which each slice should be shifted 
%                          (against the original slice position) in x- and y- direction.
%                          For each slice and for both directions (x and y) this must be given.
%                          E.g.: [0 0; 0.2 0.3; 0.5 0.5] would shift slice 0 by nothing,
%                          slice 1 by 0.2*FoV_x in x and 0.3*FoV_y in x direction and
%                          slice 3 by 0.5*FoV in both directions.
%                          size: [nSlice 2]
% 
%   Output:
% -     OutData            Shifted Output Data
% -     FoV_Phase          The Phase with which the InData gets multiplied (or in fact by exp(Phase)).
%              
%

% August 2013 - Bernhard Strasser






%% 0. Preparation

% Assign standard values to variables if nothing is passed to function.

if(~exist('InData','var'))
    display([char(10) 'You should consider inputting data which can be reconstructed. Aborting . . .'])
    return
end
if(~exist('FoV_shifts','var'))
    display([char(10) 'You have to input a shift-array if you want to shift the FoV. Aborting . . .'])
    OutData = InData;
    return
end

[nChannel nx ny nSlice nTime] = size(InData);

if(nSlice ~= size(FoV_shifts,1) || size(FoV_shifts,2) ~= 2)
    display([char(10) 'You have to input the the FoV_shifts in the proper format. Aborting . . .'])
    OutData = InData;
    return
end



%% 1. Compute the Phase to Shift the FoV


% Create a matrix with size [nChannel nx ny nSlice nTime] which linearly increases in x or y direction by -i*m*2pi, m E {-nx/2,...nx/2-1}

% Create m-Vectors
Phase_x = -nx/2:nx/2-1;
Phase_y = -ny/2:ny/2-1;

% Replicate with Slices
Phase_x = repmat(transpose(Phase_x),[1 nSlice]);
Phase_y = repmat(transpose(Phase_y),[1 nSlice]);

% Multiply with FoV Shifts and -2pi*i
for SliceIndex = 1:nSlice
    Phase_x(:,SliceIndex) = Phase_x(:,SliceIndex)*FoV_shifts(SliceIndex,1)*(-1i*2*pi);
    Phase_y(:,SliceIndex) = Phase_y(:,SliceIndex)*FoV_shifts(SliceIndex,2)*(-1i*2*pi); 
end

% Replicate to [nChannel nx ny nSlice]
Phase_x = repmat(reshape(Phase_x,[1 nx 1 nSlice]), [nChannel 1 ny 1 nTime]);
Phase_y = repmat(reshape(Phase_y,[1 1 ny nSlice]), [nChannel nx 1 1 nTime]);


FoV_Phase = Phase_x + Phase_y;





%% 2. Apply the Phase


OutData = InData .* exp(FoV_Phase);










%% ~. Postparations





