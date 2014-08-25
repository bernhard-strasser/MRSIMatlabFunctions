function [ppm_points,spec_points] = read_coord_file_0_1(file_path)
%
% read_ascconv_x_x Read ascconv header part of DICOM and Siemens raw data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function cuts out the ascconv header part of DICOM and Siemens raw data and searches for Parameters within this header. These
%
%
% [ParList,ascconv] = read_ascconv_1_2(file_path)
%
% Input: 
% -         file_path                     ...     Path of file.
%
% Output:
% -         ParList                       ...     Structure giving all the Parameters. It contains:
%           -- ParList.total_channel_no         - Number of receive-channels
%           -- ParList.ROW_raw                  - Number of measured rows (lines)
%           -- ParList.COL_raw                  - Number of measured columns (phase_encoding)
%           -- ParList.ROW                      - Number of rows (lines) AFTER zerofilling
%           -- ParList.COL                      - Number of columns (phase_encoding) AFTER zerofilling
%           -- ParList.SLC                      - Number of Slices
%           -- ParList.vecSize                  - VectorSize in spectroscopic dimension
%           -- ParList.RemoveOversampling       - Flag that determines if oversampling in spectroscopic / frequency encoding direction is removed
%           -- ParList.AsymmetricEcho           - If e.g. a GRE was measured with asymmetric echo.
%           -- ParList.InterleavedAcquisitionn  - If multiple slices were measured in interleaved acquisition mode.
%
% -         ascconv                       ...     cell array of the ascconv header: ascconv{:,1} = ParameterName, ascconv{:,2} = ParameterValue 
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





%% 5. Postparations

fclose(fid);
