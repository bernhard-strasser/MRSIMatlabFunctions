function [OutArray,ShiftMap] = FrequencyAlignment(InArray,PeakSearchSettingsOrShiftMap,ApplyAlongDim,ZerofillingFactor)
%
% FrequencyAlignment Align frequencies of csi spectra.
%
% This function was written by Bernhard Strasser, April 2015.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [OutArray,mask] = FrequencyAlignment(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
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
if(nargin < 2)
    return
end
size_InArray = size(InArray);
if(isstruct(InArray) && isfield(InArray,'mask'))
	mask = InArray.mask;
	InArray = InArray.csi;
	size_InArray = size(InArray);
	
	% Define Reference Voxel as Center of Mass Voxel
	RefVox = regionprops(mask, 'centroid');
	RefVox = round(RefVox.Centroid);
	
else
	mask = ones(size_InArray(setdiff(1:numel(size_InArray),ApplyAlongDim)));
	
	RefVox = floor(size_InArray/2)+1;
	fprintf('\n\nWARNING: No mask input for FrequencyAlignment. Reference voxel will be set as (%d,%d,%d). Might be wrong!',RefVox(1),RefVox(2),RefVox(3))
	
end



if(isfield(PeakSearchSettingsOrShiftMap,'PeakSearchPPM'))
	Settings = PeakSearchSettingsOrShiftMap;
	if(~isfield(Settings,'PolyfitRegion'))
		Settings.PolyfitRegion = [3.4 2.01];
	end
	Settings.PolyfitRegion = sort(Settings.PolyfitRegion,2,'descend');
	
else
	ShiftMap = PeakSearchSettingsOrShiftMap;
end



% 0.2 Declarations


% 0.3 Definitions
    
size_SearchArray = size_InArray; size_SearchArray(ApplyAlongDim) = size_SearchArray(ApplyAlongDim)*ZerofillingFactor;
OutArray = InArray;
SearchArray = InArray;





%% Perform Zero Order Phasing



%SearchArray = SearchArray .* exp(1i*pi);





%% 2. Exp Filtering

%SearchArray = ExponentialFilter(SearchArray, Settings.Dwelltime,12,4);




%% 1. Zerofilling & FFT SearchArray

SearchArray = Zerofilling_Spectral(SearchArray,size_SearchArray,0);
SearchArray = fftshift(fft(SearchArray,[],ApplyAlongDim),ApplyAlongDim);




%% 2. Calculate ShiftMap with Scalar Product of abs(waterpeak)

CS_vec_zf = compute_chemshift_vector_1_1(Settings.LarmorFreq,Settings.Dwelltime/10^9,Settings.vecsize*ZerofillingFactor); 
ShiftMap = zeros([size(OutArray,1) size(OutArray,2) size(OutArray,3)]);
SearchForPeak_LeftPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM));
SearchForPeak_RightPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM));
SearchForPeak_Center_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM));

ReferenceSpecMat_Spec = squeeze(SearchArray(33,33,1,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts));
ReferenceSpecMat = zeros(size(ReferenceSpecMat_Spec,1));


CircShiftVec = SearchForPeak_LeftPt_Pts-SearchForPeak_Center_Pts :1: SearchForPeak_RightPt_Pts-SearchForPeak_Center_Pts;
for i = 1:abs(SearchForPeak_RightPt_Pts-SearchForPeak_LeftPt_Pts+1)
	ReferenceSpecMat(i,:) = circshift(ReferenceSpecMat_Spec,CircShiftVec(i));
end
	
for x = 1:size(OutArray,1)
	for y = 1:size(OutArray,2)
		for z = 1:size(OutArray,3)
			
			if(mask(x,y,z) == 0 || (x==33 && y == 33 && z == 1))
				continue
			end
			
			DotProd = abs(ReferenceSpecMat) * abs(squeeze(SearchArray(x,y,z,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts)));
			ShiftMap(x,y,z) = -round( CircShiftVec(DotProd == max(DotProd)) / ZerofillingFactor); % - because we shifted the reference, but below we want to shift the other spectra
			
		end
	end
end


%% 2. Prepare searching for shift / Translate shiftmap to ShiftNoOfPointsMap

% % dwelltime gets not increased with zero_filling: zeroes get just added at the END of vector
% CS_vec_zf = compute_chemshift_vector_1_1(Settings.LarmorFreq,Settings.Dwelltime/10^9,Settings.vecsize*ZerofillingFactor); 
% 
% %SeekRef_BasisRegionWidth_SP = find(min(abs(CS_vec_zf(1) - SeekRef_BasisRegionWidth_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekRef_BasisRegionWidth_CS - CS_vec_zf));
% SearchForPeak_Center_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM));
% SearchForPeak_LeftPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM));
% SearchForPeak_RightPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM));
% PolyfitRegion_LeftPt_Pts = find(min(abs(CS_vec_zf - Settings.PolyfitRegion(1))) == abs(CS_vec_zf - Settings.PolyfitRegion(1)));
% PolyfitRegion_RightPt_Pts = find(min(abs(CS_vec_zf - Settings.PolyfitRegion(2))) == abs(CS_vec_zf - Settings.PolyfitRegion(2)));
% 
% 
% 
% %NoOfPts = SearchForPeak_RightPt_Pts - SearchForPeak_LeftPt_Pts + 1;
% %MaximumSize = size_SearchArray; MaximumSize(ApplyAlongDim) = NoOfPts;
% 
% 
% %% 3. Compute Shift
% 
% 
% SeekPeakSettings.SeekPeakRegion = [SearchForPeak_LeftPt_Pts, SearchForPeak_RightPt_Pts];
% SeekPeakSettings.PolyfitRegion = [PolyfitRegion_LeftPt_Pts, PolyfitRegion_RightPt_Pts];
% SeekPeakSettings.PolyfitOrder = 0;
% SeekPeakSettings.seek_criteria = 'findpeaks';
% SeekPeakSettings.seek_criteria_values{1} = 0.0051;
% SeekPeakSettings.seek_criteria_values{2} = 15000;
% SeekPeakSettings.seek_criteria_values{3} = 1;
% SeekPeakSettings.seek_criteria_values{4} = 3;
% SeekPeakSettings.seek_criteria_values{5} = 3;
% 
% ShiftMap = NaN([size(OutArray,1) size(OutArray,2) size(OutArray,3)]);
% 
% for x = 1:size(OutArray,1)
% 	for y = 1:size(OutArray,2)
% 		for z = 1:size(OutArray,3)
% 			
% 			if(mask(x,y,z) == 0)
% 				continue
% 			end
% 			x,y
% 			if( (x == 12 && y == 34) || (x == 41 && y == 29) || (x== 20 && y == 41) || (x == 41 && y == 45))
% 				fprintf('sdfasdf')
% 			end
% 			
% 			
% 			
% 			[NoPeakFound,Results] = seek_peak(reshape(squeeze(SearchArray(x,y,z,:)),[1 size_SearchArray(ApplyAlongDim)]),SeekPeakSettings);
% 			if(~NoPeakFound)
% 				ShiftMap(x,y,z) = Results.Peak_Point - SearchForPeak_Center_Pts;
% 			end
% 			
% 		end
% 	end
% end
% 
% 



% % Find maximum
% [Maximum,Maximum_Location] = max(real(SearchArray(:,:,:,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts)),[],ApplyAlongDim);
% % Imagine Maximum_Location = 1, i.e. the most-left point is the max. This 1 corresponds to SearchForPeak_LeftPt_Pts, i.e. 1 + SearchForPeak_LeftPt_Pts - 1 = SearchForPeak_LeftPt_Pts;
% ShiftMap = -(Maximum_Location + SearchForPeak_LeftPt_Pts - 2 - SearchForPeak_Center_Pts);  % And then we need to calculate the difference to the center, i.e. subtracting SearchForPeak_Center_Pts





%% 4. Perform Shift in Spectral Domain


OutArray = fftshift(fft(OutArray,[],ApplyAlongDim),ApplyAlongDim);

for x = 1:size(OutArray,1)
	for y = 1:size(OutArray,2)
		for z = 1:size(OutArray,3)
			
			OutArray(x,y,z,:) = circshift(squeeze(OutArray(x,y,z,:)),[ShiftMap(x,y,z) 1]);
			
		end
	end
end

OutArray = ifft(fftshift(OutArray,ApplyAlongDim),[],ApplyAlongDim);





%% 5. Postparations







