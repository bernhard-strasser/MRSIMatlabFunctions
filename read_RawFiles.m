%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%           FUNCTION TO READ IN RAW FILES          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



function image = read_RawFiles(RawFile,ROW,COL,SLC,precision)






%% 0. Preparations

Info = dir(RawFile);
if(numel(Info) < 1 || Info.bytes == 0)
	fprintf('\nScusi, but the file you gave me is zero-sized. Why should I even try to eat that?')
	image = [];
	return;
end


if(~exist('precision','var'))
	precision = 'float32';
end

if(~exist('ROW','var'))
    % determine size and compute from that ROWxROWx1
end

if(~exist('SLC','var'))
    SLC = 1;
end



%% 1. Read data

raw_fid = fopen(RawFile,'r');
image = zeros(ROW,COL,SLC);
for Slice_no = 1:SLC
	try
		image(:,:,Slice_no) = fread(raw_fid, [ROW,COL], precision);
	catch err
		image = [];
		fprintf('\nSome Error occurred when reading the raw file:\n%s',err.message)
	end
end

fclose(raw_fid);

