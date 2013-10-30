function [ParList,ascconv] = read_ascconv_1_8(file_path)
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
%           -- ParList.InterleavedSliceAcquisitionn  - If multiple slices were measured in interleaved acquisition mode.
%
% -         ascconv                       ...     cell array of the ascconv header: ascconv{:,1} = ParameterName, ascconv{:,2} = ParameterValue 
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations


% Define for which entries the ascconv should be searched for
% Search for these entries in the ascconv header part:
%if(numel(strfind(file_path, '.dat')) > 0)
    ParList_Search =  { ...
    'asCoilSelectMeas\[0\]\.asList\[\d+\]\.lRxChannelConnected',	...     % 1
    'Coil_Dumms',                                                   ...     % 27    Just a Dummy to have the parameter in the right place  
    'sRXSPEC\.alDwellTime\[\d+\]',                                  ...     % 18
    'sTXSPEC\.asNucleusInfo\[0\]\.lFrequency',                      ...     % 19	Only take the first one [0], why should there be more?
    'sKSpace\.ucDimension',                                         ...     % 8     Is this the Parameter that is different for 3D vs. 2D acquisitions????? 
    'sKSpace\.ucAsymmetricEchoAllowed',                             ...     % 11
    'sKSpace\.ucMultiSliceMode',                                    ...     % 12   
    'sKSpace\.lBaseResolution',                                     ...     % 2
    'sKSpace\.lPhaseEncodingLines',                                 ...     % 3
    'sKSpace\.lPartitions',                                         ...     % 6  
    'sSpecPara\.lFinalMatrixSizeRead',                              ...     % 4  
    'sSpecPara\.lFinalMatrixSizePhase',                             ...     % 5
    'sSpecPara\.lFinalMatrixSizeSlice',                             ...     % 34    
    'sSpecPara\.lVectorSize',                                       ...     % 9
    'sSpecPara\.ucRemoveOversampling',                              ...     % 10
    'sSpecPara.sVoI\.dPhaseFOV',                                    ...     % 23
    'sSpecPara.sVoI\.dReadoutFOV',                                  ...     % 24
    'sSpecPara.sVoI\.dThickness',                                   ...     % 25    
    'sSliceArray\.lSize',                                           ...     % 7  
    'sSliceArray.asSlice\[\d+\].dThickness',                        ...     % 20 
    'sSliceArray\.asSlice\[\d+\]\.dPhaseFOV',                       ...     % 13
    'sSliceArray\.asSlice\[\d+\]\.dReadoutFOV',                     ...     % 14
    'FOV_Partition_dumms',                                          ...     % 22	Just a Dummy to have the parameter in the right place
    'sSliceArray\.asSlice\[\d+\]\.sPosition\.dSag',                 ...     % 15    Sagittal = forehead-backhead direction = "Read"-direction (normally)
    'sSliceArray\.asSlice\[\d+\]\.sPosition\.dCor',                 ...     % 16    Coronal = left-right directions = "Phase"-direction (normally)
    'sSliceArray\.asSlice\[\d+\]\.sPosition\.dTra',                 ...     % 17    Transversal = up-down direction = "Partition or Slice"-direction (normally)
    'sSliceArray\.asSlice\[\d+\]\.sNormal\.dSag',                   ...     % 28    The x-component of the normal vector of the measured slice
    'sSliceArray\.asSlice\[\d+\]\.sNormal\.dCor',                   ...     % 29    y-component
    'sSliceArray\.asSlice\[\d+\]\.sNormal\.dTra',                   ...     % 30    z-component
    'sSliceArray.asSlice\[0].dInPlaneRot',                          ...     % 31    The InPlane (InSlice) rotation, so the rotation around the normal vector given by the upper three components
    'sGroupArray\.asGroup\[0\]\.dDistFact'                          ...     % 21
    'ucUncombImages'                                                ...     % 26 
    'sRXSPEC.lGain'                                                 ...     % 32 
    'sSpecPara\.lPhaseEncodingType'                                 ...     % 33    % 1 For Full kSpace Sampling, 2 For Elliptical Weighted Sampling, 3 for Weighted Acquisition
    };
% else
%     ParList_Search = { ...
%     'asCoilSelectMeas\[0\]\.asList\[\d+\]\.lRxChannelConnected',	...     % 1
%     'sRXSPEC.alDwellTime[0]',                                       ...     % 18
%     'sKSpace\.ucDimension',                                         ...     % 8	Is this the Parameter that is different for 3D cs. 2D acquisitions????? 
%     'sKSpace\.ucAsymmetricEchoAllowed',                             ...     % 11
%     'sKSpace\.ucMultiSliceMode',                                    ...     % 12   
%     'sKSpace\.lBaseResolution',                                     ...     % 2
%     'sKSpace\.lPhaseEncodingLines',                                 ...     % 3
%     'sKSpace\.lPartitions',                                         ...     % 6  
%     'sSpecPara\.lFinalMatrixSizeRead',                              ...     % 4  
%     'sSpecPara\.lFinalMatrixSizePhase',                             ...     % 5
%     'sSpecPara\.lVectorSize',                                       ...     % 9
%     'sSpecPara\.ucRemoveOversampling',                              ...     % 10
%     'sSliceArray\.lSize',                                           ...     % 7        
%     'sSliceArray\.asSlice\[0\]\.dPhaseFOV',                         ...     % 13
%     'sSliceArray\.asSlice\[0\]\.dReadoutFOV',                       ...     % 14
%     'sSliceArray\.asSlice\[0\]\.sPosition\.dSag',                   ...     % 15
%     'sSliceArray\.asSlice\[0\]\.sPosition\.dCor',                   ...     % 16
%     'sSliceArray\.asSlice\[0\]\.sPosition\.dTra'};                          % 17
% end


% Name the structure entries of ParList like this:
ParList_Assign = { ...
'total_channel_no_measured',                                        ...     % 1     The number of coils with which the data was measured
'total_channel_no_reco',                                            ...     % 27    The number of coils for which the data is reconstructed
'Dwelltimes',                                                       ...     % 18
'LarmorFreq',                                                       ...     % 19
'ThreeD_flag',                                                      ...     % 8
'AsymmetricEcho',                                                   ...     % 11
'InterleavedSliceAcquisition',                                      ...     % 12
'nFreqEnc',                                                         ...     % 2
'nPhasEnc',                                                         ...     % 3
'nPartEnc',                                                         ...     % 6
'nFreqEnc_FinalMatrix',                                             ...     % 4
'nPhasEnc_FinalMatrix',                                             ...     % 5
'nSLC_FinalMatrix',                                                 ...     % 34
'vecSize',                                                          ...     % 9
'RemoveOversampling',                                               ...     % 10
'VoI_Phase',                                                        ...     % 23
'VoI_Read',                                                         ...     % 24
'VoI_Partition',                                                    ...     % 25
'nSLC',                                                             ...     % 7
'SliceThickness',                                                   ...     % 20
'FoV_Phase',                                                        ...     % 13
'FoV_Read',                                                         ...     % 14
'FoV_Partition',                                                    ...     % 22
'Pos_Sag',                                                          ...     % 15
'Pos_Cor',                                                          ...     % 16
'Pos_Tra',                                                          ...     % 17
'SliceNormalVector_x',                                              ...     % 28
'SliceNormalVector_y',                                              ...     % 29
'SliceNormalVector_z',                                              ...     % 30
'InPlaneRotation',                                                  ...     % 31
'SliceGap',                                                         ...     % 21
'SaveUncombined_flag'                                               ...     % 26
'HighGain_flag'                                                     ...     % 32
'Full_ElliptWeighted_Or_Weighted_Acq'                               ...     % 33    1,2 oder 3, ob ihr wirklich...
};


% Tells function to which format it should convert the found string in the ascconv (remember: all values in the ascconv are strings):
ParList_Convert = { ...
'str2double',                                                       ...     % 1
'char',                                                             ...     % 27
'str2double',                                                       ...     % 18
'str2double',                                                       ...     % 19
'char',                                                             ...     % 8
'char',                                                             ...     % 11
'char',                                                             ...     % 12
'str2double',                                                       ...     % 2
'str2double',                                                       ...     % 3
'str2double',                                                       ...     % 6
'str2double',                                                       ...     % 4
'str2double',                                                       ...     % 5
'str2double',                                                       ...     % 34
'str2double',                                                       ...     % 9
'char',                                                             ...     % 10
'str2double',                                                       ...     % 23
'str2double',                                                       ...     % 24
'str2double',                                                       ...     % 25
'str2double',                                                       ...     % 7
'str2double',                                                       ...     % 20
'str2double',                                                       ...     % 13
'str2double',                                                       ...     % 14
'str2double',                                                       ...     % 22
'str2double',                                                       ...     % 15
'str2double',                                                       ...     % 16
'str2double',                                                       ...     % 17
'str2double',                                                       ...     % 28
'str2double',                                                       ...     % 29
'str2double',                                                       ...     % 30
'str2double',                                                       ...     % 31
'str2double',                                                       ...     % 21
'char'                                                              ...     % 26
'str2double',                                                       ...     % 32
'str2double',                                                       ...     % 33
};


% Initialize ParList
for Par_no = 1:numel(ParList_Search)
    eval([ 'ParList.' ParList_Assign{Par_no} ' = 0;' ]);
end

% open file
fid = fopen(file_path,'r');




%% 1. Track down & save element: ASCCONV

begin_found = 0;
ascconv = [];
sLine = 0;

while(sLine > -1)
    
    sLine = fgets(fid); % get string line
    
    if(not(begin_found))                                          % If begin of ascconv not yet found
        
        
        if(not(isempty(strfind(sLine,'### ASCCONV BEGIN ###'))))
            begin_found = true;                                   % If current line is begin of ascconv
        else
            continue                                              % If current line is not begin of ascconv --> read next line
        end
        
        
    else                                                          % If begin of ascconv has already been found
        
        if(not(isempty(strfind(sLine,'### ASCCONV END ###'))))    % If the end was found --> stop while loop
            break
        else
            ascconv = [ascconv; {sLine}];                         % If current line is not the end --> read in line and save it
        end
        
    end
   
        
end



%% 2. Display error & stop if no Ascconv found

if(not(begin_found))
    display(['Pfui Toifel! You gave me a file without any ascconv, I cannot digest that! Please remember that I am NOT an omnivore.' char(10) ...
             'I will stop here . . .'])
    return
end


%% 3. Convert ascconv

% Convert cell array of strings containing n x 2 entries. The first entries containing the parts before the '=' (pre-=) 
% and the second parts containing the parts after the '=' (post-=)


% Until now ascconv is a cell array of strings (with lets say 348 entries)

% This regexp makes ascconv to a cell array with 348 entries, each of these on its own a cell array of 2 strings
ascconv = regexp(ascconv, '=','split');

% This makes ascconv to a 2*348 = 696x1 cell array of strings; All odd cells contain the parts before the '=', all even cells the part after the '='
ascconv = transpose([ascconv{:}]);

% Now seperate the pre-= and the post-= parts, remove all white spaces before and after the entries.
ascconv = strtrim([ascconv(1:2:end) ascconv(2:2:end)]);

% Now we are happy and have our 348x2 cell array of strings.




%% 4. Search certain entries & Save these

% The following code performs these tasks:
% strfind(...): searches ascconv-ParameterNames (ascconv(:,1)) for the ParList_Search-strings. This results in a cell array, containing [] if in 
% the corresponding cell the Parametername was not found, and [n] if it was found in the corresponding cell on place n of the string;
% not(cellfun(...)): We then search each cell (--> cellfun) if it is empty, and negate the logical output, so that we get the non-empty cells.
% eval(...) We assign the found value to ParList.ParameterName, where ParameterName is determined by ParList_Assign. We also convert the values.

for Par_no = 1:numel(ParList_Search)
    Par_Logic = regexp(ascconv(:,1),ParList_Search{Par_no});    
    Par_Logic = not(cellfun('isempty',Par_Logic));
    if(find(Par_Logic) > 0)
        eval([ 'ParList.' ParList_Assign{Par_no} ' = ' ParList_Convert{Par_no} '(ascconv(Par_Logic,2));' ]);
    end
end


%% 5. Change & Correct certain values



% Convert from 0x1 etc. to logicals

% Remove Oversampling
if(strcmp(ParList.RemoveOversampling, '0x1')) %#ok<NODEF>
    ParList.RemoveOversampling = true;
else 
    ParList.RemoveOversampling = false;
end

% Asymmetric Echo
if(strcmp(ParList.AsymmetricEcho, '0x1'))
    ParList.AsymmetricEcho = true;
else
    ParList.AsymmetricEcho = false;
end

% Interleaved Acquisition: Slice Ordering
if(strcmp(ParList.InterleavedSliceAcquisition, '0x1'))
    ParList.InterleavedSliceAcquisition = false;
elseif(strcmp(ParList.InterleavedSliceAcquisition, '0x2'))
       ParList.InterleavedSliceAcquisition = true; 
end

% 3D_flag
if(strcmp(ParList.ThreeD_flag, '0x4'))
    ParList.ThreeD_flag = true;
else
    ParList.ThreeD_flag = false;
end

% 3D_flag
if(strcmp(ParList.SaveUncombined_flag, '0x1'))
    ParList.SaveUncombined_flag = true;
else
    ParList.SaveUncombined_flag = false;
end




% Corrections


% In the header, there is always the vector size written with oversampling removed. In the .dat-files, the oversampling is never removed. In the IMA files, it is removed, 
% if ParList.RemoveOversampling=true, otherwise not. Thus: .dat-vecsize always has to be multiplied by 2, IMA only in case of RemoveOversampling=false.
% PROBLEM: SPIRAL IS NEVER (?) OVERSAMPLED. FOR NOW: ONLY REMOVE OVERSAMPLING IF NON-3D-ACQUISITION. BUT THIS IS A HACK.
if(numel(strfind(file_path, '.dat')) > 0 || (numel(strfind(file_path, '.IMA')) > 0 && ~ParList.RemoveOversampling && ~ParList.ThreeD_flag))
 	ParList.vecSize = 2 * ParList.vecSize;   
end



% In case of Single-Slice and Multi-slice, the Partitions and the FinalMatrixSizeSlice are always set to 8, which is quite wrong.
if(~ParList.ThreeD_flag)
    ParList.nPartEnc = 1;
    ParList.nSLC_FinalMatrix = ParList.nSLC;
end



% Total channel number correction
ParList.total_channel_no_measured = numel(ParList.total_channel_no_measured);

if(numel(strfind(file_path, '.IMA')) > 0 && ~ParList.SaveUncombined_flag)
    ParList.total_channel_no_reco = 1;               % Dicom files with SaveUncombined_flag = false, are already coil combined
else
    ParList.total_channel_no_reco = ParList.total_channel_no_measured;
end

% Slice Gap
if(ParList.nSLC > 1)
    ParList.SliceGap = ParList.SliceGap .* ParList.SliceThickness(1:ParList.nSLC-1);    % There is one slice gap less than slices.
else
    ParList.SliceGap = 0;
end


% FoV in Partition Direction
ParList.FoV_Partition = sum(ParList.SliceThickness) + sum(ParList.SliceGap);





%% 5. Postparations

fclose(fid);
