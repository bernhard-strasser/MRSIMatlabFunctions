function OutArray = myrepmat(InArray,DesiredSize,ReplicationStartIndexInDesiredSize)
%
% myrepmat_1_0 Replicate matrix to desired size
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The normal MATLAB repmat function cannot add a dimension at the beginning of the array. E.g., if you have an array, size(array) = [64 64] and
% you want it to replicate it to size [32 64 64] there is no easy way to do that. If you do repmat(array, [32 1 1]) you get a array of size
% [64*32 64] = [2048 64]. This function gets you the desired [32 64 64] result.
%
%
% OutArray = myrepmat_1_0(InArray,DesiredSize,ReplicationStartIndexInDesiredSize)
%
% Input: 
% -         InArray                             ...    The input array that should be replicated
% -         DesiredSize                         ...    The desired size of the replicated array. size(InArray) must occur somewhere in this var.
% -         ReplicationStartIndexInDesiredSize  ...    Index that defines 'where size(InArray) is located in DesiredSize'. 
%                                                      E.g.: size(InArray) = [64 64], DesiredSize = [32 64 64 64 2]; ReplicationStartInd... = 2
%                                                      Tells the program, that the first two 64's are the InArray, thus InArray gets replicated
%                                                      By 32 in its 'pre-Array'-dimensions, and by [64 2] in its 'post-array'-dimensions.
%
% Output:
% -         OutArray                            ...    The replicated output array
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions



% 0.1 Preparations

OutArray = InArray;

if(~exist('DesiredSize', 'var'))
    return
end
if(~exist('ReplicationStartIndexInDesiredSize', 'var'))
    % Find first index where the size(InArray) occurs in DesiredSize and assign that value to StartIndexFor... 
    % E.g. if size(InArray) = [64 64], DesiredSize = [32 64 64 4], then we find size(InArray) in DesiredSize at starting index = 2;   
    ReplicationStartIndexInDesiredSize = strfind(reshape(DesiredSize',1,[]),size(InArray));
    if(isempty(ReplicationStartIndexInDesiredSize))
       display([ char(10) 'Error: size(InArray) must occur in DesiredSize.' ])
       return
    end
end





% 0.2 Definitions
  
% Index which determines the first 'post-array'-index
ReplicationEndIndexInDesiredSize = ReplicationStartIndexInDesiredSize + numel(size(InArray));

% In case of input-array is vector, it has only one dimension, but nevertheless matlab says it has 2 (all numbers/arrays have minimum dimensionality of 2 in MATLAB).
if(ndims(InArray) == 2 && size(InArray,2) == 1)
    ReplicationEndIndexInDesiredSize = ReplicationEndIndexInDesiredSize - 1;
end




%% 1. Replicate the 'pre-array' dimensions
% Loop over these 'pre-array' dimensions, beginning from the last one. Replicate along the dimension, so that this dimension is at the END of the
% resulting array. Then put this dimension to the beginning.
% Example: size(InArray) = [64 64]; DesiredSize = [32 2 64 64]; ReplicationStartIndexInDesiredSize = 3;
% First loop:  replicate from [64 64]    --> [64 64 2];     Reorder to [2 64 64];
% Second loop: replicate from [2 64 64]  --> [ 2 64 64 32]; Reorder to [32 2 64 64]

for Dim = ReplicationStartIndexInDesiredSize-1:-1:1                  
   
    if(DesiredSize(Dim) == 1)
        OutArray = reshape(OutArray, [1 size(OutArray)]);
    else
        Repmat_vector = [ones([1 ndims(OutArray)]) DesiredSize(Dim)];   % Create vector e.g. [1 1 1 32] for example above
        OutArray = repmat(OutArray,Repmat_vector);
        OutArray = shiftdim(OutArray,ndims(OutArray)-1);                % Shift all dimensions ndims(OutArray) positions to the left, so that the last dim gets the 1st
    end
end



%% 2. Replicate the 'post-array' dimensions
% Use repmat function directly to create all the 'post-array' dimensions, because for those the repmat works directly.
% Example: size(InArray) = [64 64]; DesiredSize = [32 64 64 1024 2]; ReplicationStartIndexInDesiredSize = 2;
% The above 'pre-array' procedure will give us: size(OutArray) = [32 64 64];
% OutArray = repmat(OutArray,[1 1 1 1024 2]). This will give us size(OutArray) = [32 64 64 1024 2]

if(ndims(OutArray) ~= numel(DesiredSize))                        % if matrix is not yet replicated to DesiredSize
    Repmat_vector = DesiredSize;
    Repmat_vector(1:ReplicationEndIndexInDesiredSize-1) = 1;

    OutArray = repmat(OutArray,Repmat_vector);
end





%% 3. Postparations

% fclose(fid)






