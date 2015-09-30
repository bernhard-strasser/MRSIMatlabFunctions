function OutArray = InterpolateCSIImages(InArray,DesiredSize,Mask,ResolutionRatio)
%
% InterpolateCSIImages Sum or Zerofill CSI Data in Image Domain
%
% This function was written by Bernhard Strasser, September 2015.
%
%
% The function sums or zerofills (in kSpace) the CSI data in image domain in order to get
% a different resolution, e.g. get the CSI from 64x64 --> 32x32.
%
%
% [OutArray,mask] = InterpolateCSIImages(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     InArray                     ...     Input array which should be interpolated
% -     DesiredSize                 ...     The desired size of the interpolated data.
%                                           You can set it to DesiredSize = 0. Then it is derived from InArray and ResolutionRatio.
% -     Mask                        ...     Mask with size(Mask) = size(InArray), so that we only sum spectra within the brain.
%                                           If Mask = 0, then set Mask = true(DesiredSize), i.e. also spectra outside of the brain are summed.
% -     ResolutionRatio             ...     The Ratio between the InArray spatial sizes and the DesiredSize. 
%                                           You can omit it or ResolutionRatio = 0. Then it is derived from InArray and DesiredSize.
%
% Output:
% -     OutArray                    ...     The interpolated output array
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 1)
	OutArray = 0;
	return
end
if(nargin < 2)
    return
end

size_InArray = size(InArray);
if(numel(DesiredSize) == 1 && DesiredSize == 0)
	DesiredSize = size_InArray(1:3) ./ ResolutionRatio;
end
if(~exist('ResolutionRatio','var') || (numel(ResolutionRatio) == 1 && ResolutionRatio == 0))
	ResolutionRatio = size_InArray(1:3) ./ DesiredSize;
end
if(~exist('Mask','var') || (numel(Mask) == 1 && Mask == 0))
	Mask = true(DesiredSize);
end
if(~islogical(Mask))
	Mask = logical(Mask);
end

% 0.2 Declarations


% 0.3 Definitions



 





%% 2. Compute Mask



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







