function InArray = hadamard_decoding(InArray,ApplyAlongDim)
%
% hadamard_decoding Perform Hadamard decoding
%
% This function was written by Gilbert Hangel, ~April 2012, and revised by Bernhard Strasser, July & August 2012.
%
%
% This function decodes hadamard-encoded data, so that you get the desired slices.
%
%
% InArray = hadamard_decoding_x(InArray,ApplyAlongDim)
%
% Input: 
% -     InArray                  ...    Input Array
% -     ApplyAlongDim            ...    Along which dimension should the hadamard decoding be applied?
%
% Output:
% -     InArray                 ...    Output Array
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 




%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
fprintf('\nHadamard decoding\t...')
tic

if(~exist('ApplyAlongDim','var'))
   ApplyAlongDim = 4; 
end


% 0.2 Declarations


% 0.3 Definitions
size_InArray = size(InArray);
  


%% 1. Memory Considerations - Find best Looping


[dummy, MemFree] = memused_linux_1_1(1);
MemNecessary = numel(InArray)*8*2*2/2^20;							% every entry of InArray is double --> 8 bytes. *2 because complex. *2 as safety measure (so InArray fits 2 times in memory,
																	% once it is already there and 2 more times it should fit in). /2^20 to get from bytes to MB.
if(numel(InArray)*8*2*2/2^20 > MemFree)
	LoopOverIndex = MemFree ./ (MemNecessary./size_InArray(1:end));	% Because the 1st index is the coil. We can never loop over the coils.
	LoopOverIndex(LoopOverIndex < 1) = NaN;
	LoopOverIndex(ApplyAlongDim) = NaN;
	LoopOverIndex = find(nanmin(LoopOverIndex));
	LoopOverIndex = LoopOverIndex(1);								% Only take the 1st if several are the minimum.
	AccessString = [	repmat(':,',[1 LoopOverIndex-1]) 'LoopIndex,' repmat(':,',[1 numel(size_InArray)-LoopOverIndex])	];
	AccessString(end) = [];
end




%% 2. Apply Hadamard Matrix

HadamardMatrix = hadamard(size_InArray(ApplyAlongDim));

if(numel(InArray)*8*2*2/2^20 > MemFree)
	
	
	for LoopIndex = 1:size_InArray(LoopOverIndex)

		tic_loop = tic;
		fprintf('\nHadamard Decoding Part %d\t...\tof %d', LoopIndex, size(InArray,LoopOverIndex))
		
		TempData = eval(['InArray(' AccessString ')']);

		% 2.1. Shift the Hadamard-Dimension to Position 1 & Reshape

		% Shift the dimensions of InArray so that the SLC-Dimension is its first. E.g.: size(InArray) = [channel ROW COL SLC vecSize] -->
		% size(ShiftedInArray) = [SLC vecSize channel ROW COL]
		TempData = shiftdim(TempData,ApplyAlongDim-1);

		% Reshape the InArray from [SLC vecSize channel ROW COL] to [SLC vecSize*channel*ROW*COL]
		ShiftedSize = size(TempData);
		TempData = reshape(TempData,[ShiftedSize(1) prod(ShiftedSize(2:end))]);

		TempData = HadamardMatrix * TempData;                       % Real Matrix Multiplication!

		% 2.3. Reshape & Shift Dimensions Back & Reshape
		% Reshape InArray from [SLC vecSize*channel*ROW*COL] to [SLC vecSize channel ROW COL]
		TempData = reshape(TempData, ShiftedSize);

		% Shift the dimensions of InArray so that size(InArray) = size_InArray. 
		TempData = shiftdim(TempData,numel(size_InArray) - (ApplyAlongDim - 1)); %#ok
		
		% To restore singleton dimensions, reshape again
		eval(['InArray(' AccessString ') = reshape(TempData,size(InArray(' AccessString ')));']);

		fprintf('\ttook\t%10.6f seconds', toc(tic_loop))       


	end

	
else
	
		% 2.1. Shift the Hadamard-Dimension to Position 1 & Reshape

		% Shift the dimensions of InArray so that the SLC-Dimension is its first. E.g.: size(InArray) = [channel ROW COL SLC vecSize] -->
		% size(ShiftedInArray) = [SLC vecSize channel ROW COL]
		InArray = shiftdim(InArray,ApplyAlongDim-1);

		% Reshape the InArray from [SLC vecSize channel ROW COL] to [SLC vecSize*channel*ROW*COL]
		ShiftedSize = size(InArray);
		InArray = reshape(InArray,[ShiftedSize(1) prod(ShiftedSize(2:end))]);

		InArray = HadamardMatrix * InArray;                       % Real Matrix Multiplication!

		% 2.3. Reshape & Shift Dimensions Back & Reshape
		% Reshape InArray from [SLC vecSize*channel*ROW*COL] to [SLC vecSize channel ROW COL]
		InArray = reshape(InArray, ShiftedSize);

		% Shift the dimensions of InArray so that size(InArray) = size_InArray. 
		InArray = shiftdim(InArray,numel(size_InArray) - (ApplyAlongDim - 1));
		InArray = reshape(InArray,size_InArray);

end


%% 3. Postparations


fprintf('\ttook\t%10.6f seconds', toc)              


















