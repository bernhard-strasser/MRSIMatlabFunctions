function [OutData,FoV_Phase] = kSpace_FoVShift(OutData, FoV_shifts) 
% 
% kSpace_FoVShift Perform In-Plane Shifts of the FoV of MR(S)I Data in kSpace
% 
%  [OutData,FoV_Phase] = kSpace_FoVShift(OutData, FoV_shifts) 
%       
%   Input:      
% -     OutData            Unshifted Data            (size: [#coils, nx, ny, nSlice, nTime]) (nTime = 1 for MRI) (For Memory reasons, same variable used for input and output data.)
%                          Must have as many slices as the OutData should have. Use similar sequence parameters as for the OutData.
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
% -     FoV_Phase          The Phase with which the OutData gets multiplied (or in fact by exp(Phase)).
%              
%

% August 2013 - Bernhard Strasser






%% 0. Preparation

% Assign standard values to variables if nothing is passed to function.

if(~exist('OutData','var'))
    display([char(10) 'You should consider inputting data which can be reconstructed. Aborting . . .'])
    return
end
if(~exist('FoV_shifts','var'))
    display([char(10) 'You have to input a shift-array if you want to shift the FoV. Aborting . . .'])
    return
end

[nChannel, nx, ny, nSlice, nTime] = size(OutData);

if(nSlice ~= size(FoV_shifts,1) || size(FoV_shifts,2) ~= 2)
    display([char(10) 'You have to input the the FoV_shifts in the proper format. Aborting . . .'])
    return
end





%% 1. Memory Considerations - Find best Looping

size_OutData = size(OutData);

[dummy, MemFree] = memused_linux(1);
MemNecessary = numel(OutData)*8*2*4/2^20;							% every entry of OutData is double --> 8 bytes. *2 because complex. *4 as safety measure (so OutData fits 2 times in memory,
																	% once it is already there and 2 more times it should fit in). /2^20 to get from bytes to MB.
ApplyAlongDims = [2 3];
if(MemNecessary > MemFree)
	LoopOverIndex = MemFree ./ (MemNecessary./size_OutData(1:end));
	LoopOverIndex(LoopOverIndex < 1) = NaN;
	LoopOverIndex(ApplyAlongDims) = NaN;
	LoopOverIndex = find(min(LoopOverIndex));
	LoopOverIndex = LoopOverIndex(1);								% Only take the 1st if several are the minimum.
	AccessString = [repmat(':,',[1 LoopOverIndex-1]) 'LoopIndex,' repmat(':,',[1 numel(size_OutData)-LoopOverIndex])];
	AccessString(end) = [];
	
	RepmatString = size_OutData;
	RepmatString(LoopOverIndex) = 1;
	RepmatString(4) = 1;
	RepmatString_x = RepmatString;
	RepmatString_x(2) = 1;
	RepmatString_y = RepmatString;
	RepmatString_y(3) = 1;
	RepmatString_x = regexprep(regexprep(num2str(RepmatString_x),' (?=\d)',','),' ','');
	RepmatString_y = regexprep(regexprep(num2str(RepmatString_y),' (?=\d)',','),' ','');
	
	
	
end







%% 2. Compute the Phase to Shift the FoV


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
% Phase_x = repmat(reshape(Phase_x,[1 nx 1 nSlice]), [nChannel 1 ny 1 nTime]);
% Phase_y = repmat(reshape(Phase_y,[1 1 ny nSlice]), [nChannel nx 1 1 nTime]);







%% 2. Repmat Apply the Phase


% Perform Noise Decorrelation by looping to avoid extensive memory usage 
if(MemNecessary > MemFree)	
	
	Phase_x = eval(['repmat(reshape(Phase_x,[1 nx 1 nSlice]),[' RepmatString_x '])']);
	Phase_y = eval(['repmat(reshape(Phase_y,[1 1 ny nSlice]),[' RepmatString_y '])']);
	FoV_Phase = Phase_x + Phase_y;
	
	for LoopIndex = 1:size(OutData,LoopOverIndex)
		TempData = eval(['OutData(' AccessString ');']);	
		TempData = TempData .* exp(FoV_Phase); %#ok
		eval(['OutData(' AccessString ') = TempData;']);
	end
	clear TempData;
	
else
	
	FoV_Phase = repmat(reshape(Phase_x,[1 nx 1 nSlice]), [nChannel 1 ny 1 nTime]) + repmat(reshape(Phase_y,[1 1 ny nSlice]), [nChannel nx 1 1 nTime]);
	OutData = OutData .* exp(FoV_Phase);

end

clear Phase_x Phase_y FoV_Phase;






%% ~. Postparations





