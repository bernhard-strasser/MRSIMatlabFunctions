function MRStruct = op_ReorderSlices(MRStruct,Settings)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%       FUNCTION TO REORDER MULTISLICE IMAGES       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Prep

if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isstruct(MRStruct))
    Bak = MRStruct; clear MRStruct; MRStruct.Data = Bak; clear Bak;
end
if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end
if(~isfield(MRStruct.RecoPar,'DataSize'))
    MRStruct.RecoPar.DataSize = size(MRStruct.Data);
end
if(~isfield(Settings,'Settings.SliceDim'))
    Settings.SliceDim = 4;
end


%% Early Exit

if(isfield_recursive(MRStruct,'RecoPar.InterleavedSliceAcquisition') && MRStruct.RecoPar.InterleavedSliceAcquisition == 0)
    return; 
end


%% Backwards compatibility:
% In general, the data can have multiple contrasts (e.g. TEs), and each could have a different image size. So best is to store those contrasts with different cell entries
% MRStruct.Data{1}, MRStruct.Data{2}, etc. However, most of my data doesn't have that structure, but MRStruct.Data is directly a matrix.
% Therefore, make it compatible with both inputs.
CellBackwardsComp_flag = false;
if(~iscell(MRStruct.Data))
    CellBackwardsComp_flag = true;
    MRStruct.Data{1} = MRStruct.Data;
    if(isfield(MRStruct,'NoiseData'))
        MRStruct.NoiseData{1} = MRStruct.NoiseData;       
    end
end


%% Reorder Data

for CurContrast = 1:numel(MRStruct.Data)

    CurData = MRStruct.Data{CurContrast};
    NoOfSlices = size(CurData, Settings.SliceDim);
    
    if(mod(NoOfSlices,2 == 0))
        SliceOrder = cat(2,2:2:NoOfSlices,1:2:NoOfSlices); % For even number of slices the slices are ordered [2 4 6 ... NoOfSlices 1 3 5 ... NoOfSlices-1]
    else
        SliceOrder = cat(2,1:2:NoOfSlices,2:2:NoOfSlices); % For odd number of slices the slices are ordered [1 3 5 ... NoOfSlices 2 4 6 ... NoOfSlices-1]
    end
    [dum,SliceOrderedInd] = sort(SliceOrder);
    
    
    % This is a cryptic way of doing: MRStruct.Data{CurContrast} = CurData(:,:,...,SliceOrderedInd,:,:,...); where the dimension-position of SliceOrderedInd is exactly
    % Settings.SliceDim. But since the Slice dimension can change, we cannot just write e.g. Out = In(:,:,:,SliceOrderedInd,:,:,:); What if the slice dimension is not
    % dimension 4 but 5? Below code is general enough to cope with that.
    SubStruct.type = '()'; SubStruct.subs = repmat({':'},1,ndims(CurData));
    SubStruct.subs{Settings.SliceDim} = SliceOrderedInd;
    MRStruct.Data{CurContrast} = subsref(CurData,SubStruct);
    if(isfield(MRStruct,'NoiseData'))
        MRStruct.NoiseData{CurContrast} = subsref(MRStruct.NoiseData{CurContrast},SubStruct);        
    end

end
    

%% Backwards compatibility:
% In general, the data can have multiple contrasts (e.g. TEs), and each could have a different image size. So best is to store those contrasts with different cell entries
% MRStruct.Data{1}, MRStruct.Data{2}, etc. However, most of my data doesn't have that structure, but MRStruct.Data is directly a matrix.
% Therefore, make it compatible with both inputs.
if(CellBackwardsComp_flag)
    MRStruct.Data = MRStruct.Data{1};
    if(isfield(MRStruct,'NoiseData'))
        MRStruct.NoiseData = MRStruct.NoiseData{1};       
    end
end



%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);




