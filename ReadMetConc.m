function [MetConc, MetNames, Extra, ExtraNames] = ReadMetConc(csv_file,table_file)
% ReadMetConc reads the metabolite concentrations of a CSV-file created by LCModel
% Usage: 
% Input: Path of CSV-file; 
% Output: Metabolite Concentrations (array), Metabolite Names (array)
%
%
%
% ##################################################
% ########### Function to Hell and Back ############
% ##################################################
%            Written by Bernhard Strasser




%% 0. Preps

if(~exist('csv_file','var'))
	csv_file = 0;
end
if(~exist('table_file','var'))
	table_file = 0;
end


%% 1. Read CSV-file

if(~isempty(csv_file) && ischar(csv_file))
	MetConc_dummy = importdata(csv_file);
	if(~isempty(MetConc_dummy))
		MetConc = MetConc_dummy.data(3:end);                % The first two columns in file are called "row" and "col" for processing MRSI data at once, don't need that info here
		MetNames = MetConc_dummy.textdata(3:end);           % The first two columns in file are called "row" and "col" for processing MRSI data at once, don't need that info here
	else
		MetConc = [];
		MetNames = [];
	end
end







%% 2. Read table file

if(~isempty(table_file) && ischar(table_file))
	tablefid = fopen(table_file,'r');
	table = textscan(tablefid,'%s');
	table = table{:};
	SNR = table(find(logical(cellfun(@numel,regexpi(table, 'S/N')))) + 2);
	LW = table(find(logical(cellfun(@numel,regexpi(table, 'FWHM')))) + 2);
	fclose(tablefid);
	if(isempty(SNR))
		SNR{1} = '0';
	end
	if(isempty(LW))
		LW{1} = '0';
	end
	Extra = cat(2,str2num(SNR{1}),str2num(LW{1}));
	ExtraNames = {'SNR','LW'};
	
end







