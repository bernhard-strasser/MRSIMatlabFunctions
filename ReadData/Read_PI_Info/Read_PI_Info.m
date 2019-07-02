function R_Total = Read_PI_Info(path)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, Dec 2014.
%
%
% The function reads in the header info of a jmrui txt file.
%
%
% Data = read_JmruiTxtData(file)
%
% Input: 
% -         file                    ...     jmrui txt file.
%
% Output:
% -         Data                    ...     The read in header info.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 





%% 0. Preparation

if(~exist('path','var'))
	fprintf('\nError: No file passed over.\n')
	R_Total = 0;
	return
end
if(~exist(path,'file'))
	fprintf('\nError: File\n%s\ndoes not exist.\n',path)
	R_Total = 0;
	return
end


%% 1. Gather information from header




fiddy = fopen(path,'r');


tline = fgetl(fiddy);

while ischar(tline)
	
	isRline = regexp(tline,'R\t\t\t\t\t=\t','ONCE');
	if(~isempty(isRline))
		R_Total = regexp(tline,'\t=\t.*','match');
		R_Total = regexprep(R_Total,'\t=\t','');
		R_Total = str2double(R_Total{:});
	end
	tline = fgetl(fiddy);

	
end


%% 


fclose(fiddy);





