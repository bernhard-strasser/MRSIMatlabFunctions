%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%           FUNCTION TO READ IN RAW FILES          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



function write_RawFiles(RawArray,RawFile,precision)






%% 0. Preparations


if(~exist('precision','var'))
	precision = 'float32';
end

Directory = regexp(RawFile,'/'); Directory = RawFile(1:Directory(end)-1);
if(~exist(Directory,'dir'))
	mkdir(Directory);
end




%% 1. Write data

raw_fid = fopen(RawFile,'w');
for Slice_no = 1:size(RawArray,3)
    fwrite(raw_fid, RawArray(:,:,Slice_no), precision);
end
fclose(raw_fid);

