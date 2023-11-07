function OutData = ShiftSpectrum(OutData, ShiftMap_Pts,ApplyAlongDim) 
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
if(~exist('ShiftMap_Pts','var'))
    display([char(10) 'You have to input a shift-array if you want to shift the FoV. Aborting . . .'])
    return
end


sizenumelwas4 = false;
if(numel(size(OutData)) == 4)
	sizenumelwas4 = true;
	OutData = reshape(OutData,[1 size(OutData)]); ApplyAlongDim = ApplyAlongDim + 1;
end
[nChannel, nx, ny, nSlice, nTime] = size(OutData);
size_OutData = size(OutData);
size_Comp = size_OutData(2:end-1); size_Comp(size_Comp == 1) = [];

if(sum(size_Comp ~= size(ShiftMap_Pts)) > 0 )
    display([char(10) 'You have to input the the FoV_shifts in the proper format. Aborting . . .'])
    return
end



	






%% 1. Memory Considerations - Find best Looping


size_OutData = size(OutData);

[dummy, MemFree] = memused_linux(1);
MemNecessary = numel(OutData)*8*3*2/2^20;						% every entry of OutData is double --> 8 bytes. *2 because complex. *3 as safety measure (so OutData fits 2 times in memory,
PossibleLoopIndices = [1];										% once it is already there and 2 more times it should fit in). /2^20 to get from bytes to MB.
if(MemNecessary > MemFree)
	LoopOverIndex = MemFree ./ (MemNecessary./size_OutData(PossibleLoopIndices));
	LoopOverIndex(LoopOverIndex < 1) = NaN;
	LoopOverIndex = PossibleLoopIndices(LoopOverIndex == nanmin(LoopOverIndex));
	LoopOverIndex = LoopOverIndex(1);							% Only take the 1st if several are the minimum. Add 1 because the channel index (first) is not considered
	AccessString = [repmat(':,',[1 LoopOverIndex-1]) 'LoopIndex,' repmat(':,',[1 numel(size_OutData)-LoopOverIndex])];
	AccessString(end) = [];
end






%% 2. Compute the Phase to Shift the spectrum


% Create a matrix with size [nChannel nx ny nSlice nTime] which linearly increases in t direction by -i*m*2pi, m E {-nx/2,...nx/2-1}

% Repmat ShiftMap_Pts
ShiftMap_Pts = repmat(ShiftMap_Pts,[1 1 1 nTime]);

% Create Phase Vector
Phase_t = (0:nTime-1)*(-1i*2*pi/nTime);

% Replicate 
Phase_t = repmat(reshape(Phase_t,[1 1 1 nTime]),[nx ny nSlice 1]);

% Multiply with FoV Shifts and -2pi*i
Phase_t = Phase_t.*ShiftMap_Pts; clear ShiftMap_Pts;






%% 2. Repmat & Apply the Phase


% Perform Noise Decorrelation by looping to avoid extensive memory usage 
if(MemNecessary > MemFree)	
		
	for LoopIndex = 1:size(OutData,LoopOverIndex)
		TempData = eval(['OutData(' AccessString ');']);	
		TempData = TempData .* exp(Phase_t); %#ok
		eval(['OutData(' AccessString ') = TempData;']);
	end
	clear TempData;
	

	
	
	
	
	
else
	Phase_t = repmat(reshape(Phase_t,[1 size(Phase_t)]), [nChannel 1 1 1 1]);
	OutData = OutData .* exp(Phase_t);

end

clear Phase_t Phase_y FoV_Phase;






%% ~. Postparations

if(sizenumelwas4)
	OutData = reshape(OutData,size_OutData(2:end));
end



