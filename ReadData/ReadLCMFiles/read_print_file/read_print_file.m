function PlotData = read_print_file(file_path,BeginStr,EndStr)
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


%% 1. 


Text = fileread(file_path);

blocksstart = regexp(Text,BeginStr);
Text = Text(blocksstart(1):end);

blockend = regexp(Text,EndStr);
Text = Text(1:blockend(1));


filebyline = regexp(Text, '\n', 'split');
filebyline = filebyline(2:end);
%remove empty lines
filebyline( cellfun(@isempty,filebyline) ) = [];

% Remove all lines at end without numbers, until we find a line with numbers
for CurLine = numel(filebyline):-1:1
    if(isempty(regexp(filebyline{CurLine},'\d\.\d')))   % Doesnt contain a decimal point
        filebyline(CurLine) = [];
    else
        break
    end
end

% concat lines if no char is found in beginning
LinesToCat = cellfun(@(v)v(2:5),filebyline,'UniformOutput',false);
LinesToCat = cellfun(@isempty,regexp(LinesToCat,'[a-zA-Z]'));
LinesToCat = sort(find(LinesToCat),'descend');

for catline = LinesToCat
    if(catline == 1)
        continue
    end
    filebyline{catline-1} = strcat(2,filebyline{catline-1},' ', filebyline{catline});
    filebyline(catline) = [];
    fprintf('')
end


%split by fields
filebyfield = regexp(filebyline, '\s+', 'split');
%pad out so each line has the same number of fields
numfields = cellfun(@length, filebyfield);
maxfields = max(numfields);
fieldpattern = repmat({[]}, 1, maxfields);
firstN = @(S,N) S(1:N);
filebyfield = cellfun(@(S) firstN([S,fieldpattern], maxfields), filebyfield, 'Uniform', 0);
%switch from cell vector of cell vectors into a 2D cell
fieldarray = vertcat(filebyfield{:});

% LinesWithNumbers = find(~isnan(cellfun(@str2double,fieldarray(:,2))));
% LinesWithNumbers = sort(LinesWithNumbers,'descend');
% 
% for catline = transpose(LinesWithNumbers)
%    fieldarray(catline-1,:) = cat(2,fieldarray(catline-1,:),fieldarray(catline,:));
%    fieldarray(catline,:) = [];
%    fprintf('')
% end


%% 5. Postparations

asdf
