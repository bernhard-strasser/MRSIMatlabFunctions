function PlotData = read_coord_file(file_path)
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



% open file
fid = fopen(file_path,'r');
ppm_found = false;
sLine = 0;


SearchSet = {'NY phased data points follow','NY points of the fit to the data follow','NY background values follow'};       % Search for these strings
IdentifyWithSet = {'Spec','Fit','Baseline'};                                                                                % And if we find one, identify with this data set
EndOfReadSearchSet = {'\d+ lines in following diagnostic table:'};                                                          % If we find one of these strings, stop reading


%% 1. Read ppm-scale


while(sLine > -1)
    
    sLine = fgets(fid); % get string line

    if(not(isempty(strfind(sLine,'points on ppm-axis '))))
        
        ppm_found = true;                                               % If current line is begin of ppm-axis
        PlotData.NumberOfPoints = textscan(sLine,'%d');                 % Find out the number of points included in the ppm-scale which is written in this line.
        PlotData.NumberOfPoints = double(PlotData.NumberOfPoints{:});   % Convert from cell to double.
        
        % Read ppm-points
        PlotData.ppm_points = fscanf(fid,'%f',PlotData.NumberOfPoints);
        
        break

    else
        continue                                              % If current line is not begin of ppm-axis --> read next line
    end
        
end



% Display error & stop if no ppm scale found

if(not(ppm_found))
    fprintf(['\nPfui Toifel! You gave me a wrong file, I cannot digest this file:\n%s\n' ...
    'Please remember that I am NOT an omnivore. I will stop here . . .\n'],file_path)
    PlotData.ppm_points = 0;
    fclose(fid);
	return
end




%% Read rest


while(sLine > -1)
    
    sLine = fgets(fid); % get string line
    
    % stop reading if we find one of the strings in EndOfReadSearchSet
    if(sum(~cellfun(@isempty,regexp(sLine,EndOfReadSearchSet))) > 0)
        break
    end
    
    % The data set defining lines need to include alphabetic characters and do not have numbers, except for 'Conc. = +xyz.abcd'
    ConcIndex = regexp(sLine,'Conc\. = .?\d+\.\d+');
    ConcFound = numel(ConcIndex) > 0;
    if(numel(regexp(sLine,'[a-zA-Z]+')) == 0 || (numel(regexp(sLine,'\d+')) > 0 && ~ConcFound))
        continue
    end    
    % Otherwise it seems to be a line defining a data set 
    
    
    % See if we can search one of the phrases defined in SearchSet
    CurDataSet = sLine;
    if(~ConcFound)
        SearchSet_FoundIndex = ~cellfun(@isempty,regexp(CurDataSet,SearchSet));

        % If we found something, replace with IdentifyWithSet
        if(sum(SearchSet_FoundIndex) > 0)
            CurDataSet = IdentifyWithSet{SearchSet_FoundIndex};
        else                                                                % It's neither a 'Metabo Conc. = ...', nor one of the phrases in SearchSet
            error('\nError in read_coord_file: I found a line which does not fit in my patterns. Abort.')   
        end
    else
        CurDataSet = CurDataSet(1:ConcIndex(1)-1);              % Remove all this concentration stuff
    end
    
    % In any case, remove all whitespace characters, and check if string is non-empty
    CurDataSet = regexprep(CurDataSet,'\s','');
    % Replace invalid field names
    CurDataSet = regexprep(CurDataSet,'2HG','TwoHG');
    if(~isempty(regexp(CurDataSet(1),'[^a-zA-Z]','ONCE')))
        CurDataSet(1) = [];
    end
   
    % Read data
    if(~isempty(CurDataSet))
        if(ConcFound)
            PlotData.Metabos.(CurDataSet) = fscanf(fid,'%f',PlotData.NumberOfPoints);
        else
            PlotData.(CurDataSet) = fscanf(fid,'%f',PlotData.NumberOfPoints);            
        end
    end
    
    

    
end






%% 5. Postparations

fclose(fid);
