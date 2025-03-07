function OutArray = op_XFillOrCutData(OutArray,Settings)
%
% ZerofillOrCutkSpace Zerofill or Cut k-Space Data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function either zerofills or cuts data in k-space. For zerofilling time data 
% (zerofilling only one-sided instead of on both sides like in k-space), use "Zerofilling_Spectral".
%
%
% OutArray = ZerofillOrCutkSpace(OutArray,Settings.Zerofill_To,Settings.PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied.
% -     Settings.Zerofill_To         ...    Array to which data should be zerofilled or cut. E.g. size(OutArray) = [32 64 64 512], Settings.Zerofill_To = [32 128 128 512]. 
% -     Settings.PerformFFT_flag     ...    If it is 1, the image gets Fourier transformed to k-space before applying the filter, 
%                                           and transformed back to image domain afterwards
% -     Settings.X                   ...    What value should be zerofilled. E.g. X=0 will zero-fill, X=1, will one-fill
% -     Settings.AppendZerosTo       ...    For each dimension specify if data should be appended to beginning, end or both, by inputting a cell. E.g.
%                                           {'Beginning','End','Both','Both'} appends zeros at beginning in 1st dim, end on 2nd dim, and both sides for 3rd and 4th
%                                           Useful e.g. for zerofilling Partial-Fourier data.
%
% Output:
% -     OutArray                     ...    The filtered/masked output array
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
NoStruct_flag = false;
if(~isstruct(OutArray))
    NoStruct_flag = true;
    Bak = OutArray;
    clear OutArray;
    OutArray.Data = Bak;
    clear Bak;
end
if(~isfield(Settings,'PerformFFT_flag'))
   Settings.PerformFFT_flag = false; 
end 
if(~isfield(Settings,'X'))
    Settings.X = 0;      % Default: Zerofilling
end

if(~isfield(Settings,'AppendZerosTo'))
    Settings.AppendZerosTo = repmat({'Both'},[1 numel(Settings.Zerofill_To)]);
end
if(numel(Settings.AppendZerosTo) < numel(Settings.Zerofill_To))
    Settings.AppendZerosTo = cat(2,Settings.AppendZerosTo,repmat({'Both'},[1 numel(Settings.Zerofill_To)-numel(Settings.AppendZerosTo)]));
end

% 0.2 Declarations


% 0.3 Definitions
    
size_OutArray = size(OutArray.Data);
size_OutArray = [size_OutArray ones([1 numel(Settings.Zerofill_To)-numel(size_OutArray)])];
AppendZeros = round(Settings.Zerofill_To - size_OutArray);
%AppendZeros(AppendZeros < 0) = 0;
ApplyAlongDims = find(ne(AppendZeros,0));



 





%% 1. FFT to k-space

if(Settings.PerformFFT_flag)

    for filter_dim = ApplyAlongDims
        OutArray.Data = ifftshift(OutArray.Data,filter_dim);
        OutArray.Data = ifft(OutArray.Data,[],filter_dim);
        OutArray.Data = fftshift(OutArray.Data,filter_dim);
    end  

end




%% 2. Compute Mask


for dummy_dim = ApplyAlongDims

	if(AppendZeros(dummy_dim) > 0)
        
        if(strcmpi(Settings.AppendZerosTo{dummy_dim},'Beginning'))
            NoOfXsAtBeginning = AppendZeros(dummy_dim);
            NoOfXsAtEnd = 0;
        elseif(strcmpi(Settings.AppendZerosTo{dummy_dim},'End'))
            NoOfXsAtBeginning = 0;
            NoOfXsAtEnd = AppendZeros(dummy_dim);            
        else
            NoOfXsAtBeginning = ceil(AppendZeros(dummy_dim)/2);
            NoOfXsAtEnd = floor(AppendZeros(dummy_dim)/2);            
        end
        
		AppendZeros_AtBeginning = myrepmat(Settings.X,[size_OutArray(1:dummy_dim-1) NoOfXsAtBeginning size_OutArray(dummy_dim+1:end)]);
		AppendZeros_AtEnd = myrepmat(Settings.X,[size_OutArray(1:dummy_dim-1) NoOfXsAtEnd size_OutArray(dummy_dim+1:end)]);
		OutArray.Data = cat(dummy_dim,AppendZeros_AtBeginning,OutArray.Data,AppendZeros_AtEnd);
	else
		center = floor(size(OutArray.Data,dummy_dim)/2) + 1;
		AccessPart = {num2str(center - floor(Settings.Zerofill_To(dummy_dim)/2)), num2str(center + ceil(Settings.Zerofill_To(dummy_dim)/2) - 1)};
		AccessString = [repmat(':,',[1 dummy_dim-1]) AccessPart{1} ':' AccessPart{2} ',' repmat(':,',[1 numel(size_OutArray)-dummy_dim])];
		AccessString(end) = [];
		OutArray.Data = eval(['OutArray.Data(' AccessString ');']);
	end
    size_OutArray = size(OutArray.Data);
    
end


%% 4. FFT to ImageSpace


if(Settings.PerformFFT_flag)
    
    for filter_dim = ApplyAlongDims
        OutArray.Data = ifftshift(OutArray.Data,filter_dim);
        OutArray.Data = fft(OutArray.Data,[],filter_dim);
        OutArray.Data = fftshift(OutArray.Data,filter_dim);
    end  
    
end




%% 5. Postparations


if(NoStruct_flag)
   Bak = OutArray.Data;
   clear OutArray;
   OutArray = Bak;
end




