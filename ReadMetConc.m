function [MetConc, MetNames] = ReadMetConc(csv_file)
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



%% 1. Read CSV-file

MetConc_dummy = importdata(csv_file);

if(~isempty(MetConc_dummy))
    MetConc = MetConc_dummy.data(3:end);                % The first two columns in file are called "row" and "col" for processing MRSI data at once, don't need that info here
    MetNames = MetConc_dummy.textdata(3:end);           % The first two columns in file are called "row" and "col" for processing MRSI data at once, don't need that info here
else
    MetConc = [];
    MetNames = [];
end


   
