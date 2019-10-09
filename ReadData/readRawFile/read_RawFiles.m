function image = read_RawFiles(RawFile,ROW,COL,SLC,precision)
% read_RawFiles Read a simple binary 'raw' file, e.g. created by mnc2raw
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function 
%
%
% image = read_RawFiles(RawFile,ROW,COL,SLC,precision)
%
% Input: 
% -         RawFile                       ...     Path of file.
% -         ROW                           ...     The size of the image in the first direction.
% -         COL                           ...     The size of the image in the second direction.
% -         SLC                           ...     The size of the image in the third direction.
% -         precision                     ...     The precision of the doubles in the file given as a string. E.g. if you used
%                                                 minctoraw bla.mnc -nonormalize -float > bla.raw, precision is 'float32'
%                                                 Further examples: ???
%                                                 Default: 'float32'
%
% OR:
% -         RawFile                       ...     Path of file.
% -         [ROW,COL,SLC]                 ...     The size of the image in the first direction.
% -         precision                     ...     The precision of the doubles in the file given as a string. E.g. if you used
%                                                 minctoraw bla.mnc -nonormalize -float > bla.raw, precision is 'float32'
%                                                 Further examples: ???
%                                                 Default: 'float32'%
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

if(numel(ROW) == 3)
   if(exist('COL','var'))
       precision = COL;
   end
   SLC = ROW(3);
   COL = ROW(2);
   ROW = ROW(1);
   
end

Info = dir(RawFile);
if(numel(Info) < 1 || Info.bytes == 0)
	fprintf('\nError in read_RawFiles: Scusi, but the file you gave me is zero-sized. Why should I even try to eat that?\n')
	image = [];
	return;
end


if(~exist('precision','var'))
	precision = regexpi(RawFile,'_prec_');
	if(isempty(precision))
		precision = 'float32';
	else
		precision = regexp(RawFile(precision:end),'_prec_[a-zA-Z]+[0-9]*','match');
		precision = precision{1}(7:end);

	end
end

if(~exist('ROW','var'))
	ROW = regexpi(RawFile,'_x');
	if(isempty(ROW))
		ROW = '128';
	else
		ROW = regexp(RawFile(ROW:end),'_x[0-9]+','match'); 
		ROW = str2double(ROW{1}(3:end));
	end
end

if(~exist('COL','var'))
	COL = regexpi(RawFile,'_y');
	if(isempty(COL))
		COL = '128';
	else
		COL = regexp(RawFile(COL:end),'_y[0-9]+','match'); 
		COL = str2double(COL{1}(3:end));
	end
end

if(~exist('SLC','var'))
	SLC = regexpi(RawFile,'_z');
	if(isempty(SLC))
		SLC = '128';
	else
		SLC = regexp(RawFile(SLC:end),'_z[0-9]+','match'); 
		SLC = str2double(SLC{1}(3:end));
	end
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

