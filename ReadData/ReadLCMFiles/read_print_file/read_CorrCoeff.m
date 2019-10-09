function Output = read_CorrCoeff(file_path)
%
% read_coord_file Read Coordinate LCModel file
%
% This function was written by Bernhard Strasser, September 2017.
%
%
% The function reads "Coordinate" LCModel files. In these files the plots of the output-ps-files of LCModel are written in numbers, e.g. to create your own plots. 
%
%
% PlotData = read_coord_file(file_path)
%
% Input: 
% -         file_path                     ...     Path of file.
%
% Output:
% -         PlotData                      ...     Structure containing all the plots. Example fields:
%                   .NumberOfPoints       ...     Number of points that were read in.
%                   .ppm_points           ...     Read in points of the ppm-axis
%                   .Spec                 ...     Read in points of Input spectrum
%                   .Fit                  ...     Read in points of LCModel fit
%                   .Baseline             ...     Read in points of the fitted baseline
%                   .Metabos              ...     Structure containing all the metabolite fits as subfields:
%                           .xyz          ...     Read in points of metabolite xyz
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations

if(~exist('EndStr','var'))
    EndStr = '\n\n\n';
end
    
% % open file
% fid = fopen(file_path,'r');
% ppm_found = false;
% sLine = 0;
% 
% 
% SearchSet = {'NY phased data points follow','NY points of the fit to the data follow','NY background values follow'};       % Search for these strings
% IdentifyWithSet = {'Spec','Fit','Baseline'};                                                                                % And if we find one, identify with this data set
% EndOfReadSearchSet = {'\d+ lines in following diagnostic table:'};                                                          % If we find one of these strings, stop reading


%% 1. Read file, split in lines, remove lines, concat lines

Text = read_print_file(file_path,'Correlation coefficients',[],true);
if(isempty(Text))
    Output = [];
    return;
end

% Split lines
filebyline = regexp(Text, '\n', 'split');

% Remove first and empty lines
filebyline = filebyline(2:end);
filebyline( cellfun('isempty',filebyline) ) = [];

% Remove last lines which is just repetition of first lines
% At the end the metabolites are listed again, but we have this info already at the top.
% Therefore, remove all lines at end without numbers, until we find a line with numbers
for CurLine = numel(filebyline):-1:1
    if(isempty(regexp(filebyline{CurLine},'\d\.\d','ONCE')))   % Doesnt contain a decimal point
        filebyline(CurLine) = [];
    else
        break
    end
end


% Concat lines
% Unfortunately, long lines are split into several lines. Need to undo this.
% Since we know that all normal lines start with the metabolite name in first column, flag those lines that don't start with some characters in front as being
% extension lines, that should be concatenated to the previous one.
% concat lines if no char is found in beginning
LinesToCat = cellfun(@(v)v(2:5),filebyline,'UniformOutput',false);
LinesToCat = cellfun('isempty',regexp(LinesToCat,'[a-zA-Z]'));
LinesToCat = sort(find(LinesToCat),'descend');
for catline = LinesToCat
    if(catline == 1)
        continue
    end
    filebyline{catline-1} = strcat(2,filebyline{catline-1},' ', filebyline{catline});
    filebyline(catline) = [];
end


%% Replace strange characters, Remove Heading and Trailing Whitespace

filebyline = regexprep(filebyline,char(2),' ');
filebyline = cellfun(@strtrim,filebyline,'UniformOutput',false);

filebyline{1} = [' ',filebyline{1}];    % The first line with all the metabos need to be shifted one to right


%% Split in columns, remove columns

%split by fields
filebyfield = regexp(filebyline, '\s+', 'split');
%pad out so each line has the same number of fields
numfields = cellfun(@length, filebyfield);
maxfields = max(numfields);
fieldpattern = repmat({[]}, 1, maxfields);
firstN = @(S,N) S(1:N);
filebyfield = cellfun(@(S) firstN([S,fieldpattern], maxfields), filebyfield, 'Uniform', 0);
%switch from cell vector of cell vectors into a 2D cell
Output = vertcat(filebyfield{:});

% Remove empty columns
EmptyCells = cellfun('isempty',Output);
EmptyCells(~EmptyCells) = cellfun('isempty',regexp(Output(~EmptyCells),['[^\s^' char(2) ']']));
EmptyCells = all(EmptyCells);
Output(:,EmptyCells) = [];


%% Replace empty entries, convert strings to double

% Replace empty entries with '0.0'-strings
EmptyCells = cellfun('isempty',Output);
EmptyCells(1,1) = false;    % The first one is ok if it is emtpy
Output(EmptyCells) = {'0.0'};

% Convert string numbers to double
Output(2:end,2:end) = cellfun(@(x){str2double(x)},(Output(2:end,2:end)));


%% 5. Postparations


