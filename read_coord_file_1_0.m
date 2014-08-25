function [ppm_points,spec_points] = read_coord_file_1_0(file_path)
%
% read_coord_file_x_y Read Coordinate LCModel file
%
% This function was written by Bernhard Strasser, September 2012.
%
%
% The function reads "Coordinate" LCModel files. In these files the plots of the output-ps-files of LCModel are written in numbers, e.g. to create your own plots. 
%
%
% [ppm_points,spec_points] = read_coord_file_x_y(file_path)
%
% Input: 
% -         file_path                     ...     Path of file.
%
% Output:
% -         ppm_points                    ...     The chemical shift scale of the spectrum.
% -         spec_points                   ...     These are the phased spectrum points.
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
spec_found = false;
ppm_points = 0;
spec_points = 0;
sLine = 0;




%% 1. Read ppm-scale


while(sLine > -1)
    
    sLine = fgets(fid); % get string line

    if(not(isempty(strfind(sLine,'points on ppm-axis '))))
        
        ppm_found = true;                                     % If current line is begin of ppm-axis
        NumberOfPoints = textscan(sLine,'%d');                % Find out the number of points included in the ppm-scale which is written in this line.
        NumberOfPoints = double(NumberOfPoints{:});           % Convert from cell to double.
        
        % Read ppm-points
        ppm_points = fscanf(fid,'%f',NumberOfPoints);
        
        break

    else
        continue                                              % If current line is not begin of ppm-axis --> read next line
    end
        
end



% Display error & stop if no ppm scale found

if(not(ppm_found))
    display(['Pfui Toifel! You gave me a wrong file, I cannot digest that! Please remember that I am NOT an omnivore.' char(10) ...
             'I will stop here . . .'])
    return
end






%% 2. Read Spectrum

if(ppm_found)
    

    while(sLine > -1)

        sLine = fgets(fid); % get string line

        if(not(isempty(strfind(sLine,'NY phased data points follow'))))   % If current line is begin of spectrum
            
            spec_found = true;                                   
            
            % Read spec-points
            spec_points = fscanf(fid,'%f',NumberOfPoints);

            break

        else
            continue                                                      % If current line is not begin of ascconv --> read next line
        end

    end   
    

end



% Display error & stop if no spectrum points found

if(not(spec_found))
    display(['Pfui Toifel! You gave me a wrong file, I cannot digest that! Please remember that I am NOT an omnivore.' char(10) ...
             'I will stop here . . .'])
    return
end





%% 3. Postparations

fclose(fid);
