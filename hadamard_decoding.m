function OutArray = hadamard_decoding_9(InArray,ApplyAlongDim)
%
% hadamard_decoding Perform Hadamard decoding
%
% This function was written by Gilbert Hangel, ~April 2012, and revised by Bernhard Strasser, July & August 2012.
%
%
% This function decodes hadamard-encoded data, so that you get the desired slices.
%
%
% OutArray = hadamard_decoding_x(InArray,ApplyAlongDim)
%
% Input: 
% -     InArray                  ...    Input Array
% -     ApplyAlongDim            ...    Along which dimension should the hadamard decoding be applied?
%
% Output:
% -     OutArray                 ...    Output Array
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
    
OutArray = InArray;






%% 1. Shift the Hadamard-Dimension to Position 1 & Reshape

% Shift the dimensions of InArray so that the SLC-Dimension is its first. E.g.: size(InArray) = [channel ROW COL SLC vecSize] -->
% size(ShiftedInArray) = [SLC vecSize channel ROW COL]
OutArray = shiftdim(OutArray,ApplyAlongDim-1);

% Reshape the OutArray from [SLC vecSize channel ROW COL] to [SLC vecSize*channel*ROW*COL]
ShiftedOutArray_Size = size(OutArray);
OutArray = reshape(OutArray,[ShiftedOutArray_Size(1) prod(ShiftedOutArray_Size(2:end))]);





%% 2. Apply Hadamard Matrix

HadamardMatrix = hadamard(size(InArray,ApplyAlongDim));
OutArray = HadamardMatrix * OutArray;                       % Real Matrix Multiplication!






%% 4. Reshape & Shift Dimensions Back & Reshape

% Reshape OutArray from [SLC vecSize*channel*ROW*COL] to [SLC vecSize channel ROW COL]
OutArray = reshape(OutArray, ShiftedOutArray_Size);

% Shift the dimensions of OutArray so that size(OutArray) = size(InArray). 
OutArray = shiftdim(OutArray,numel(size(InArray)) - (ApplyAlongDim - 1));

% To restore singleton dimensions, reshape again
OutArray = reshape(OutArray, size(InArray));






%% 5. Postparations

