function HeaderInfo = read_JmruiTxtHeader(file)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, Dec 2014.
%
%
% The function reads in the header info of a jmrui txt file.
%
%
% HeaderInfo = read_JmruiTxtHeader(file)
%
% Input: 
% -         file                    ...     jmrui txt file.
%
% Output:
% -         HeaderInfo              ...     The read in header info.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 





%% 0. Preparation

if(~exist('file','var'))
	fprintf('\nError: No file passed over.\n')
	HeaderInfo = 0;
	return
end



%% 1. Gather information from header





fid = fopen(file,'r');
tline = fgetl(fid);

if(~strcmpi(tline,'jMRUI Data Textfile'))
	fprintf('\nError: File does not seem to be jMRUI txt file.\n')
	fclose(fid);
	HeaderInfo = 0;
	return
end	

repet = 0;
while (~strcmp(tline, 'Signal and FFT'))
	repet = repet + 1;
	if(repet > 50)
		fprintf('\nError: File does not seems to be corrupt.\n')
		fclose(fid);
		HeaderInfo = 0;
		return
	end
		
	if(numel(strfind(tline, 'PointsInDataset:'))>0)
		HeaderInfo.total_points = str2double(strrep(tline, 'PointsInDataset: ', ''));
	elseif(numel(strfind(tline, 'SamplingInterval:'))>0)
		HeaderInfo.dwelltime = str2double(strrep(tline, 'SamplingInterval: ', ''))/1000;   %/1000 so that the unit is s
	elseif(numel(strfind(tline, 'TransmitterFrequency:'))>0)
		HeaderInfo.Frequency = str2double(strrep(tline, 'TransmitterFrequency: ', ''))/1000000;   % [MHz]
	end
	tline = fgetl(fid);
end
fclose(fid);



