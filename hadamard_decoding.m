function OutArray = hadamard_decoding(InArray,ApplyAlongDim)
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
	
	Mini = min([LoopOverIndex,ApplyAlongDim]);
	Diffi = abs(LoopOverIndex - ApplyAlongDim);
	First = {'i','LoopIndex'}; Second = First{(Mini == LoopOverIndex) + 1}; First(strcmpi(First,Second)) = []; First = First{:};
	AccessString = [	repmat(':,',[1 Mini-1]) First repmat(':,',[1 Diffi]) Second repmat(':,',[1 numel(size_InArray)-(Mini-1)-(Diffi-1)-2])	];
	AccessString(end) = [];
	
else
	AccessString = [repmat(':,',[1 ApplyAlongDim-1]) 'i,' repmat(':,',[1 numel(size_InArray)-ApplyAlongDim])];
	AccessString(end) = [];
end

AccessString_j = regexprep(AccessString,'i','j');





%% 2. Apply Hadamard Matrix

HadamardMatrix = hadamard(size(InArray,ApplyAlongDim));
OutArray = zeros(size_InArray);

for i = 1:size(HadamardMatrix,1)
	TempData = 0;
	for j = 1:size(HadamardMatrix,2)
		
		if(numel(InArray)*8*2*2/2^20 > MemFree)
			for LoopIndex = 1:size_InData(LoopOverIndex)
				TempData = TempData + HadamardMatrix(i,j) * eval(['InArray(' AccessString_j ')']);
			end
		else
			TempData = TempData + HadamardMatrix(i,j) * eval(['InArray(' AccessString_j ')']);
		end
		
	end
	eval(['OutArray(' AccessString ') = TempData;']);
end






%% 3. Postparations

