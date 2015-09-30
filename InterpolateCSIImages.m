function OutArray = InterpolateCSIImages(InArray,Mask,DesiredSize,ResolutionRatio)
%
% InterpolateCSIImages Sum or Zerofill CSI Data in Image Domain
%
% This function was written by Bernhard Strasser, September 2015.
%
%
% The function sums or zerofills (in kSpace) the CSI data in image domain in order to get
% a different resolution, e.g. get the CSI from 64x64 --> 32x32
%
%
% [OutArray,mask] = InterpolateCSIImages(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                           Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                           Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                           and transformed back to image domain afterwards
%
% Output:
% -     OutArray                     ...     The filtered/masked output array
% -     mask                         ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 1)
	OutArray = 0;
	return
end
if(nargin < 3)
    return
end
size_InArray = size(InArray);
if(numel(DesiredSize) == 1 && DesiredSize == 0)
	DesiredSize = size_InArray(1:3) ./ ResolutionRatio;
end
if(~exist('ResolutionRatio','var') || (numel(ResolutionRatio) == 1 && ResolutionRatio == 0))
	ResolutionRatio = size_InArray(1:3) ./ DesiredSize;
end
if(~islogical(Mask))
	Mask = logical(Mask);
end

% 0.2 Declarations


% 0.3 Definitions



 









%% 2. Compute Mask

% OutArray = InArray;
% for dim = 1:numel(ResolutionRatio)
% 		
% 		
% 		if(ResolutionRatio(dim) == 1)
% 			continue;
% 		end
% 		
% 		InArray = zeros([size_InArray(1:dim-1) size_InArray(dim)/ResolutionRatio(dim) size_InArray(dim+1:end)]);
% 		if(ResolutionRatio(dim) > 1)
% 			% Sum
% 			for XorYorZ_new = 1:size(InArray,dim)
% 				XorYorZ_old_dwn = ResolutionRatio(dim)*(XorYorZ_new - 1) + 1;
% 				XorYorZ_old_up = XorYorZ_new*ResolutionRatio(dim);
% 				AssignStr = [ 'InArray(' repmat(':,',[1 dim-1]) num2str(XorYorZ_new) repmat(',:',[1 numel(size_InArray)-dim]) ') = ' ];
% 				AssignStr = [ AssignStr 'sum(OutArray(' repmat(':,',[1 dim-1]) num2str(XorYorZ_old_dwn) ':' num2str(XorYorZ_old_up) repmat(',:',[1 numel(size_InArray)-dim]) '),' num2str(dim) ');' ];
% 				eval(AssignStr);
% 			end
% 			InArray = InArray / ResolutionRatio(dim);
% 			
% 		elseif(ResolutionRatio(dim) < 1)
% 			% Zerofill in kSpace
% 		end
% 		OutArray = InArray;
% 		size_InArray = size(InArray);
% end



OutArray = zeros([DesiredSize prod(size_InArray(4:end))]);
for z_new = 1:DesiredSize(3)
	for y_new = 1:DesiredSize(2)
		for x_new = 1:DesiredSize(1)
			
			x_old_dwn = ResolutionRatio(1)*(x_new - 1) + 1;
			x_old_up = ResolutionRatio(1)*x_new;
			y_old_dwn = ResolutionRatio(2)*(y_new - 1) + 1;
			y_old_up = ResolutionRatio(2)*y_new;
			z_old_dwn = ResolutionRatio(3)*(z_new - 1) + 1;
			z_old_up = ResolutionRatio(3)*z_new;
			
			% The mask which should be potentially summed
			CurMask = false(size_InArray(1:3));
			CurMask(x_old_dwn:x_old_up,y_old_dwn:y_old_up,z_old_dwn:z_old_up) = true;
			CurMask = logical(CurMask .* Mask);
			NoOfSummedVoxels = sum(sum(sum(CurMask)));
			
			if(NoOfSummedVoxels <= 0)
				continue
			end
			
			CurMask = myrepmat(CurMask,size_InArray);
			OutArray(x_new,y_new,z_new,:) = sum(reshape(InArray(CurMask),[NoOfSummedVoxels prod(size_InArray(4:end))])) / NoOfSummedVoxels;
			
		end
	end
end

OutArray = reshape(OutArray,[DesiredSize size_InArray(4:end)]);







%% 5. Postparations







