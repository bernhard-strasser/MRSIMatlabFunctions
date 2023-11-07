function write_RawFiles(RawArray,RawFile,precision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%           FUNCTION TO WRITE IN RAW FILES          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% write_RawFiles(RawArray,RawFile,precision)




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
fwrite(raw_fid, RawArray, precision);


fclose(raw_fid);

