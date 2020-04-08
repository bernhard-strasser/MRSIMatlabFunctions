function image = ReadMncFile(MncFile)
% read_MncFiles Read a simple binary 'raw' file, e.g. created by mnc2raw
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function 
%
%
% image = read_MncFiles(MncFile,ROW,COL,SLC,precision)
%
% Input: 
% -         MncFile                       ...     Path of file.
%
%
% Output:
% -        image                          ...     The data which was read in.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None








%% 0. Preparations

Info = dir(MncFile);
if(numel(Info) < 1 || Info.bytes == 0)
	fprintf('\nError in ReadMncFile: Scusi, but the file you gave me is zero-sized. Why should I even try to eat that?\n')
	image = [];
	return;
end

% Unfortunately, this doesn't work. Minc needs to be setup before starting matlab! Only then it works!
% if(exist('MncPath','var'))
%     % If it is already loaded, you can still use this function, even if the following command gives an error
%     [tmp, tmpout] = unix(['. ' MncPath]);    
% end


%% 1. Find Data Size and Data Dimension Order

[tmp,tmpout] = unix(['mincinfo ' MncFile]);

if(tmp ~= 0)
    fprintf('\nWarning from ReadMncFile: Could not run mincinfo on file "%s".\nMaybe minc-tools is not installed?',MncFile)
    image = [];
    return
end

xspace = regexp(tmpout,'xspace\s*\d+','match'); xspace = str2double(regexp(xspace{1},'\d+','match'));
yspace = regexp(tmpout,'yspace\s*\d+','match'); yspace = str2double(regexp(yspace{1},'\d+','match'));
zspace = regexp(tmpout,'zspace\s*\d+','match'); zspace = str2double(regexp(zspace{1},'\d+','match'));

% Somehow it seems that the order written in a mnc file is always [z y x], even if in the header the first dimension appears as x...
% xspace_Ind = regexp(tmpout,'xspace\s*\d+');
% yspace_Ind = regexp(tmpout,'yspace\s*\d+');
% zspace_Ind = regexp(tmpout,'zspace\s*\d+');
% [~,DimOrder] = sort([xspace_Ind, yspace_Ind, zspace_Ind]);  % DimOrder = [3 2 1] e.g. means x is 3rd dim, y is second, and z is 1st
DimOrder = [3 2 1];
DimSizes = [xspace, yspace, zspace];
DimSizes = DimSizes(DimOrder);


%% 2. Convert File to Raw

% Generate random unique id:
temp =  java.util.UUID.randomUUID;
myuuid = char(temp.toString);

RawFile = regexprep(MncFile,'.mnc',['_' myuuid '.raw']);
[tmp,tmpout] = unix(['minctoraw -nonormalize -float ' MncFile ' > ' RawFile]);



%% 3. Read Raw File

raw_fid = fopen(RawFile,'r');

image = fread(raw_fid, 'float');
image = reshape(image,DimSizes);

% Reorder image to be always x,y,z, even if mnc was e.g. z,y,x
image = permute(image,DimOrder);

fclose(raw_fid);


%% 4. Delete Raw File

delete(RawFile);



