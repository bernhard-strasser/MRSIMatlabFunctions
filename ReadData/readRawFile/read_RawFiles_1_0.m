%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%           FUNCTION TO READ IN RAW FILES          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



function image = read_RawFiles_1_0(RawFile,precision,ROW,COL,SLC)






%% 0. Preparations


pause on


if(~exist('ROW','var'))
    % determine size and compute from that ROWxROWx1
end

if(~exist('SLC','var'))
    SLC = 1;
end



%% 1. Read data

raw_fid = fopen(RawFile,'r');
image = fread(raw_fid, [ROW,COL], precision);
fclose(raw_fid);

