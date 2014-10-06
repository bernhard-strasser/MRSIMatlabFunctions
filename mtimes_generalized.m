function C = mtimes_generalized(A,B, CommonDims)
%
% mtimes_generalized Generalize matrixmultiplication to general Arrays. 
%
% This function was written by Bernhard Strasser, August 2014.
%
%
% The function takes an Array 
% A of size a1 x a2 x ... x an x k      and 
% B of size k x b1 x b2 x ... x bm
% and computes Array C of size a1 x a2 x ... x an x b1 x b2 x ... x bm       by 
%
%
% C = mtimes_generalized(A,B, CommonDims)
%
% Input: 
% -     A           ...    The left side Array of A*B. size(A,end) = size(B,1) must hold!
% -     B           ...    The right side Array of A*B. size(A,end) = size(B,1) must hold!
% -     CommonDims  ...    If A and B have a common dimension which should not be replicated for B,
%                          e.g. if size(A) = [2 3 4], size(B) = [4 2 5]. Normally this would lead to
%                          size(C) = [2 3 2 5], but if you state that the 2 is a common dimension
%                          by CommonDims = [1     2], it will be replicated to size(C) = [2 3 5].
%                                           ^     ^
%                                           |     |
%                    place of common dim in A  |  place of common dim in B.
%                          CommonDims can be a whole list with size(CommonDims) = [n 2].
%                    
%
% Output:
% -     C           ...    C = A*B (array multiplication, NOT element-wise multiplication).
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 2)
	error('Err:NotEnoughInput','Not enough input arguments.')
end
if(size(A,numel(size(A))) ~= size(B,1))
	error('Err:ArraySizeIncons','Error using mtimes_generalized:\nInner array dimensions must agree\n(size(A,end) = size(B,1)).')
end
if(~exist('CommonDims','var') || sum(sum(CommonDims <= 0)) > 0)
	CommonDims = [];
end
if( numel(CommonDims) > 0 &&   (size(CommonDims,2) < 2 || size(CommonDims,1) > numel(size(A)) || size(CommonDims,1) > numel(size(B)))   )
	warning('Warn:InconCommonDims','Warning: Your CommonDims variable is inconsistent (with A and B)!\n')
	CommonDims = [];
end	
	

% 0.2 Declarations


% 0.3 Definitions
sizeA = size(A);
sizeB = size(B);


%% 1. Figure out the DimensionCorrespondence

if(numel(CommonDims) > 1)
	sizeB(CommonDims(:,2)) = [];
	AddDims = 1:numel(sizeB)+size(CommonDims(:,2)); AddDims = setdiff(AddDims, CommonDims(:,2));
	Corry = [ zeros([1 numel(sizeA)-1])   AddDims ];
	Corry(CommonDims(:,1)) = CommonDims(:,2); 
else
	Corry = [ zeros([1 numel(sizeA)-1])   1:numel(sizeB) ];	
end




%% 2. Replicate, multiply & sum

A = myrepmat(A,[sizeA(1:end-1) sizeB], [1:numel(sizeA) zeros([1 numel(sizeB)-1])]);
B = myrepmat(B,[sizeA(1:end-1) sizeB], Corry);

C = A; clear A;
C = C .* B;
C = sum(C,numel(sizeA));
C = squeeze_single_dim(C,numel(sizeA));





%% 5. Postparations







