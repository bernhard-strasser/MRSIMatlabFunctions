function Text = read_print_file(file_path,BeginStr,EndStr,quiet_flag)
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

if(~exist('EndStr','var') || isempty(EndStr))
    EndStr = '\n\n\n';
end
if(~exist('quiet_flag','var'))
    quiet_flag = false;
end
    
%
if(~exist(file_path,'file'))
    if(~quiet_flag)
        fprintf('\nError in read_print_file: File %s does not exist. Abort.',file_path)
    end
    Text = [];
    return;
end

%% 1. 


Text = fileread(file_path);


if(exist('BeginStr','var') && ~isempty(BeginStr))
    blocksstart = regexp(Text,BeginStr);
    if(isempty(blocksstart))
        if(~quiet_flag)
            fprintf('\nError in read_print_file: File %s does not contain string %s. Abort.',file_path,BeginStr)
        end
        Text = [];
        return;
    end
    Text = Text(blocksstart(1):end);
end
    
if(exist('EndStr','var') && ~isempty(EndStr))
    blockend = regexp(Text,EndStr);
    Text = Text(1:blockend(1));
end





%% 5. Postparations


