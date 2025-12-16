function [ParList,ascconv] = read_ascconv_VE11_eh(file_path,NumberOfAscconvEnds)
%
% read_ascconv Read ascconv header part of DICOM and Siemens raw data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function cuts out the ascconv header part of DICOM and Siemens raw data and searches for Parameters within this header. These
%
%
% [ParList,ascconv] = read_ascconv(file_path)
%
% Input: 
% -         file_path                     ...     Path of file.
%
% Output:
% -         ParList                       ...     Structure giving all the Parameters. It contains among many others:
%
%           -- ParList.total_channel_no_measured         - Number of receive-channels that were measured            
%           -- ParList.total_channel_no_reco             - Number of receive-channels that were are in the file (DICOM)             
%           -- ParList.Dwelltimes                        - Dwelltime             
%           -- ParList.LarmorFreq                        - LarmorFrequency                
%           -- ParList.ThreeD_flag                       - flag telling you if measurement is real 3D measurement or 2D/Multislice                     
%           -- ParList.AsymmetricEcho                    - If AssymetricEcho was allowed                    
%           -- ParList.InterleavedSliceAcquisition       - If Slices were not measured consecutively like 1,2,3,4,... but e.g. 1,3,2,4         
%           -- ParList.nFreqEnc                          - Number of measured points in frequency encoding direction       
%           -- ParList.nPhasEnc                          - Number of measured points in phase encoding direction           
%           -- ParList.nPartEnc                          - Number of measured points in partition encoding direction (only for 3D)
%           ...
%           ...
%           ...
%
%
% -         ascconv                       ...     cell array of the ascconv header: ascconv{:,1} = ParameterName, ascconv{:,2} = ParameterValue 
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations

if(~exist('NumberOfAscconvEnds','var'))
    NumberOfAscconvEnds = 0;
end


% Define for which entries the ascconv should be searched for
% Search for these entries in the ascconv header part:
ParList_Search =  { ...
'sCoilSelectMeas\.aRxCoilSelectData\[0\]\.asList\[\d+\]\.lRxChannelConnected',	...     % 1
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
'sSpecPara\.sVoI\.sPosition\.dSag',								...     % 43    Sagittal = forehead-backhead direction = "Read"-direction (normally)
'sSpecPara\.sVoI\.sPosition\.dCor',								...     % 44    Coronal = left-right directions = "Phase"-direction (normally)
'sSpecPara\.sVoI\.sPosition\.dTra',								...     % 45    Transversal = up-down direction = "Partition or Slice"-direction (normally)
'sSpecPara\.sVoI\.sNormal\.dSag',								...     % 46    The x-component of the normal vector of the measured slice
'sSpecPara\.sVoI\.sNormal\.dCor',								...     % 47    y-component
'sSpecPara\.sVoI\.sNormal\.dTra',								...     % 48    z-component
'sSpecPara\.sVoI\.dInPlaneRot',									...     % 49    The InPlane (InSlice) rotation, so the rotation around the normal vector given by the upper three components
'sSliceArray\.asSlice\[\d+\]\.sPosition\.dSag',                 ...     % 15    Sagittal = forehead-backhead direction = "Read"-direction (normally)
'sSliceArray\.asSlice\[\d+\]\.sPosition\.dCor',                 ...     % 16    Coronal = left-right directions = "Phase"-direction (normally)
'sSliceArray\.asSlice\[\d+\]\.sPosition\.dTra',                 ...     % 17    Transversal = up-down direction = "Partition or Slice"-direction (normally)
'sSliceArray\.asSlice\[\d+\]\.sNormal\.dSag',                   ...     % 28    The x-component of the normal vector of the measured slice
'sSliceArray\.asSlice\[\d+\]\.sNormal\.dCor',                   ...     % 29    y-component
'sSliceArray\.asSlice\[\d+\]\.sNormal\.dTra',                   ...     % 30    z-component
'sSliceArray.asSlice\[0].dInPlaneRot',                          ...     % 31    The InPlane (InSlice) rotation, so the rotation around the normal vector given by the upper three components
'sGroupArray\.asGroup\[0\]\.dDistFact',                         ...     % 21
'ucUncombImages',                                               ...     % 26 
'sRXSPEC.lGain',                                                ...     % 32 
'sSpecPara\.lPhaseEncodingType',                                ...     % 33    % 1 For Full kSpace Sampling, 2 For Elliptical Weighted Sampling, 3 for Weighted Acquisition
'alTE\[\d+\]',                                                  ...     % 35
'sWiPMemBlock\.tFree',														...		% 36
'(?<!a)lAverages',												...		% 37	Only search for lAverages excluding alAverages, because this belongs to 'sDiffusion.alAverages.__attribute__.size'
'sWipMemBlock\.alFree\[(\d){1,2}\]',							...		% 38	All variables set in Special Card + those from above
'sKSpace.ucPhasePartialFourier',								...		% 39
'sKSpace.ucSlicePartialFourier',								...		% 40
'tProtocolName',												...		% 41	How the "Sequence" at the scanner when it was measured was named
'tSequenceFileName',											...		% 42	Which sequence was used
'alTI\[\d+\]',                                                  ...		% 51	TI
'alTR\[\d+\]',                                                  ...		% 52	TR
'sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude',            ...     % 50
'lTotalScanTimeSec',                                            ...     % 53    Total scan time
'sTXSPEC\.asNucleusInfo\[0\]\.tNucleus',                        ...     % 54    Nucleus
'sAdjData\.sAdjVolume.dThickness',                              ...     % 55    AdjBoxSize slice
'sAdjData\.sAdjVolume.dPhaseFOV',                               ...     % 56    AdjBoxSize Phase
'sAdjData\.sAdjVolume.dReadoutFOV',                             ...     % 57    AdjBoxSize Read
'sWipMemBlock\.adFree\[(\d){1,2}\]'                     		...		% 58	All variables set in Special Card + those from above
};



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
'PosVOI_Sag',                                                       ...     % 43
'PosVOI_Cor',                                                       ...     % 44
'PosVOI_Tra',                                                       ...     % 45
'SliceNormalVector_VOI_x',                                          ...     % 46
'SliceNormalVector_VOI_y',                                          ...     % 47
'SliceNormalVector_VOI_z',                                          ...     % 48
'InPlaneRotation_VOI',                                              ...     % 49
'Pos_Sag',                                                          ...     % 15
'Pos_Cor',                                                          ...     % 16
'Pos_Tra',                                                          ...     % 17
'SliceNormalVector_x',                                              ...     % 28
'SliceNormalVector_y',                                              ...     % 29
'SliceNormalVector_z',                                              ...     % 30
'InPlaneRotation',                                                  ...     % 31
'SliceGap',                                                         ...     % 21
'SaveUncombined_flag',                                              ...     % 26
'HighGain_flag',                                                    ...     % 32
'Full_ElliptWeighted_Or_Weighted_Acq',                              ...     % 33    1,2 oder 3, ob ihr wirklich...
'TEs',                                                              ...     % 35
'wipMemBlock_tFree',												...		% 36
'nAverages',														...		% 37
'WipMemBlock_alFree',												...		% 38
'PhasePartialFourier',												...		% 39
'SlicePartialFourier',												...		% 40
'tProtocolName',													...		% 41
'tSequenceFileName',												...		% 42
'TIs',                                                              ...		% 51
'TR',                                                               ...		% 52
'RefAmplitude',                                                     ...     % 50
'TotalScanTime',                                                    ...     % 53
'Nucleus',                                                          ...     % 54
'AdjBoxSize_Slice',                                                 ...     % 55
'AdjBoxSize_Phase',                                                 ...     % 56
'AdjBoxSize_Read',                                                  ...     % 57
'WipMemBlock_adFree'												...		% 58
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
'str2double',                                                       ...     % 43
'str2double',                                                       ...     % 44
'str2double',                                                       ...     % 45
'str2double',                                                       ...     % 46
'str2double',                                                       ...     % 47
'str2double',                                                       ...     % 48
'str2double',                                                       ...     % 49
'str2double',                                                       ...     % 15
'str2double',                                                       ...     % 16
'str2double',                                                       ...     % 17
'str2double',                                                       ...     % 28
'str2double',                                                       ...     % 29
'str2double',                                                       ...     % 30
'str2double',                                                       ...     % 31
'str2double',                                                       ...     % 21
'char',                                                             ...     % 26
'str2double',                                                       ...     % 32
'str2double',                                                       ...     % 33
'str2double',                                                       ...     % 35
'char',																...		% 36
'str2double',														...		% 37
'str2double',														...		% 38
'char',																...		% 39
'char',																...		% 40
'char',																...		% 41
'char',																...		% 42
'str2double',                                                       ...     % 51
'str2double',                                                       ...     % 52
'str2double',														...		% 50
'str2double',														...		% 53
'char',                                                             ...     % 54
'str2double',                                                       ...     % 55
'str2double',														...		% 56
'str2double',														...		% 57
'str2double'														...		% 58
};



BasicInfo = io_GetBasicTwixfileInfos(file_path);
fid = fopen(file_path,'r','n','UTF-8');
wiff = fread(fid,BasicInfo.BeginningOfData,'int8=>char');   % Read in only until beginning of data. Otherwise the search of big files takes ages.
wiff = convertCharsToStrings(wiff);
Tmp = regexp(wiff,'<ParamDouble."flReadoutOSFactor">  { <Precision> 6  \d\.\d*\s*}','match');
if(isempty(Tmp))        % For new XA versions
    Tmp = regexp(wiff,'<ParamDouble."flReadoutOSFactor">{\s*\d*\s*}','match','once');
end
% [~, Tmp] = unix(['grep -a -m 4 ReadoutOSFactor ' file_path]);
% Tmp = regexp(Tmp, '<ParamDouble."flReadoutOSFactor">  { <Precision> 6  .*  }','match');
Tmp = regexp(Tmp{1}, '{.*}','match');
Tmp = regexprep(Tmp,'<Precision> 6  ','');
Tmp = regexprep(Tmp,'{\s+','');
Tmp = regexprep(Tmp,'\s+}','');
ReadoutOSFactor = str2double(Tmp{1});
if(isempty(ReadoutOSFactor) || isnan(ReadoutOSFactor))
    ReadoutOSFactor = 2;
end



% Initialize ParList
for Par_no = 1:numel(ParList_Search)
    eval([ 'ParList.' ParList_Assign{Par_no} ' = ' ParList_Convert{Par_no} '(''0'');' ]);
end
ParListBak = ParList;

% rewind
frewind(fid);


firstInt  = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');

NScans = secondInt;
measID = fread(fid,1,'uint32');
fileID = fread(fid,1,'uint32');
measOffset = cell(1, NScans);
measLength = cell(1, NScans);
for k=1:NScans
    measOffset{k} = fread(fid,1,'uint64');
    measLength{k} = fread(fid,1,'uint64'); 
    fseek(fid, 152 - 16, 'cof');
end


%% 1. Track down & save element: ASCCONV


for CurSet = 1:NScans
    
    ParList = ParListBak;   % For now only take last ParList. If others are also necessary, can make a cell out of it later
    CurAscconvEndCount = 1;
    MeaningfulAscconvNotFound = true;
    fseek(fid,measOffset{CurSet},'bof');
    while(MeaningfulAscconvNotFound)

begin_found = 0;
ascconv = [];
sLine = 0;
Ascconv_No = 0;

while(sLine > -1)
    
    sLine = fgets(fid); % get string line

    CurPos = ftell(fid);
    if(CurPos > BasicInfo.BeginningOfData)
        break;
    end

    if(not(begin_found))                                            % If begin of ascconv not yet found
        
        
        if(not(isempty(strfind(sLine,'### ASCCONV BEGIN'))))
            begin_found = true;                                     % If current line is begin of ascconv
            Ascconv_No = Ascconv_No+1;
        else
            continue                                                % If current line is not begin of ascconv --> read next line
        end
        
        
    else                                                            % If begin of ascconv has already been found
        
        if(not(isempty(strfind(sLine,'### ASCCONV END ###'))))      % If the end was found --> stop while loop
            if(Ascconv_No < NumberOfAscconvEnds) 		% THIS IS A HACK!
                begin_found = false;
                ascconv = [];
            else
                break
            end
        elseif(not(isempty(strfind(sLine,'### ASCCONV BEGIN'))))    % In very rare cases there are two 'ASCCONV BEGIN's in the header...
            ascconv = [];
        else
            ascconv = [ascconv; {sLine}];                           % If current line is not the end --> read in line and save it
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
if(mod(numel(ascconv),2) == 1)
    ascconv = ascconv(1:end-1); 
end
ascconv = strtrim([ascconv(1:2:end) ascconv(2:2:end)]);

% Now we are happy and have our 348x2 cell array of strings.







%% 4. Search certain entries & Save these

% The following code performs these tasks:
% strfind(...): searches ascconv-ParameterNames (ascconv(:,1)) for the ParList_Search-strings. This results in a cell array, containing [] if in 
% the corresponding cell the Parametername was not found, and [n] if it was found in the corresponding cell on place n of the string;
% not(cellfun(...)): We then search each cell (--> cellfun) if it is empty, and negate the logical output, so that we get the non-empty cells.
% eval(...) We assign the found value to ParList.ParameterName, where ParameterName is determined by ParList_Assign. We also convert the values.

for Par_no = 1:numel(ParList_Search)
    TokenAssignString = regexp(ascconv(:,1),ParList_Search{Par_no},'tokens');    
    Par_Logic = not(cellfun('isempty',TokenAssignString));
	TokenAssignString = TokenAssignString(Par_Logic); 
	while(iscell([TokenAssignString{:}]))
		TokenAssignString = [TokenAssignString{:}];
	end

	if(numel(TokenAssignString) > 0)
        TokenAssignString = cell2mat(strcat(TokenAssignString,'+1,'));
		TokenAssignString = ['([' TokenAssignString(1:end-1) '])'];
	else
        TokenAssignString = '';
	end
	if(find(Par_Logic) > 0)
        eval([ 'ParList.' ParList_Assign{Par_no} TokenAssignString ' = ' ParList_Convert{Par_no} '(ascconv(Par_Logic,2));' ]);
	end
end



        %%
        if((ParList.TR == 0 || all(ParList.TEs == 0) || all(ParList.Dwelltimes == 0) || strcmpi(ParList.Nucleus,'0')) && CurAscconvEndCount < 10)
            CurAscconvEndCount = CurAscconvEndCount + 1;
            continue
        else
            MeaningfulAscconvNotFound = false;
        end
    end


    
    
    %% 5. Change & Correct certain values


% FoV in Partition Direction
if(ParList.nPartEnc == 1)
    ParList.FoV_Partition(1) = ParList.VoI_Partition(1);
else
    ParList.FoV_Partition = sum(ParList.SliceThickness) + (ParList.nSLC > 1)*sum(ParList.SliceGap);
end

% Convert from 0x1 etc. to logicals

% Remove Oversampling
if(~isempty(regexp(ParList.RemoveOversampling, '1(?!.)','once')))			% (?!.) means not followed by anything (anything = .). E.g. 0x10 will not be matched, but 0x1 will.
    ParList.RemoveOversampling = true;
else 
    ParList.RemoveOversampling = false;
end

% Asymmetric Echo
if(~isempty(regexp(ParList.AsymmetricEcho, '1(?!.)','once')))
    ParList.AsymmetricEcho = true;
else
    ParList.AsymmetricEcho = false;
end

% Interleaved Acquisition: Slice Ordering
if(~isempty(regexp(ParList.InterleavedSliceAcquisition, '1(?!.)','once')))
	ParList.InterleavedSliceAcquisition = false;
elseif(~isempty(regexp(ParList.InterleavedSliceAcquisition, '2(?!.)','once')))
	ParList.InterleavedSliceAcquisition = true;
elseif(~isempty(regexp(ParList.InterleavedSliceAcquisition, '4(?!.)','once')))
	ParList.InterleavedSliceAcquisition = 3;				% Single-Shot	   
end

% 3D_flag
if(~isempty(regexp(ParList.ThreeD_flag, '4(?!.)','once')))
    ParList.ThreeD_flag = true;
else
    ParList.ThreeD_flag = false;
end

% 3D_flag
if(~isempty(regexp(ParList.SaveUncombined_flag, '1(?!.)','once')))
    ParList.SaveUncombined_flag = true;
else
    ParList.SaveUncombined_flag = false;
end

% Phase Partial Fourier
if(~isempty(regexp(ParList.PhasePartialFourier, '1(?!.)','once')))
    ParList.PhasePartialFourier = 4/8;
elseif(~isempty(regexp(ParList.PhasePartialFourier, '2(?!.)','once')))
    ParList.PhasePartialFourier = 5/8;
elseif(~isempty(regexp(ParList.PhasePartialFourier, '4(?!.)','once')))
    ParList.PhasePartialFourier = 6/8;
elseif(~isempty(regexp(ParList.PhasePartialFourier, '8(?!.)','once')))
    ParList.PhasePartialFourier = 7/8;
else
    ParList.PhasePartialFourier = 0;
end

% Slice Partial Fourier
if(~isempty(regexp(ParList.SlicePartialFourier, '1(?!.)','once')))
    ParList.SlicePartialFourier = 4/8;
elseif(~isempty(regexp(ParList.SlicePartialFourier, '2(?!.)','once')))
    ParList.SlicePartialFourier = 5/8;
elseif(~isempty(regexp(ParList.SlicePartialFourier, '4(?!.)','once')))
    ParList.SlicePartialFourier = 6/8;
elseif(~isempty(regexp(ParList.SlicePartialFourier, '8(?!.)','once')))
    ParList.SlicePartialFourier = 7/8;
else
    ParList.SlicePartialFourier = 0;
end

% Check which sequence we are assuming. This is implemented really badly by checking the SequenceFileName, SequenceDescription, and SequenceString.
% Better would be a unique code for each sequence
ParList = CheckAssumedSequence(ParList,file_path);


    
% In the header, there is always the vector size written with oversampling removed. In the .dat-files, the oversampling is never removed. In the IMA files, it is removed, 
% if ParList.RemoveOversampling=true, otherwise not. Thus: .dat-vecsize always has to be multiplied by 2, IMA only in case of RemoveOversampling=false.
% PROBLEM: SPIRAL IS NEVER (?) OVERSAMPLED. FOR NOW: ONLY REMOVE OVERSAMPLING FOR FULLY SAMPLED DATA SET. BUT THIS IS A HACK.

ParList.ReadoutOSFactor = ReadoutOSFactor;

if(numel(strfind(file_path, '.dat')) > 0 ...
	|| (numel(strfind(file_path, '.IMA')) > 0 && ~ParList.RemoveOversampling && ParList.Full_ElliptWeighted_Or_Weighted_Acq ~= 1))
 	ParList.vecSize = ParList.ReadoutOSFactor * ParList.vecSize;
    if( strcmp(ParList.AssumedSequence,'ViennaCRT') || strcmp(ParList.AssumedSequence,'AntoinesEccentricOrRosette') )      % if Rollercoaster
        ParList.Dwelltimes = ParList.ReadoutOSFactor * ParList.Dwelltimes;    % Because we have oversampling enabled, and the scanner interprets that we oversample in the spectral domain, which we dont...
    end
end
if( numel(strfind(file_path, '.IMA') > 0) ) 
	if( strcmp(ParList.AssumedSequence,'BorjanSpiral'))	% PREVIOUSLY HAD: ParList.Full_ElliptWeighted_Or_Weighted_Acq ~= 4 &&. BUT SOMETIMES CSI IS ALSO WEIGHTED!!! % WITH THAT I ASSUME THAT 
	 	ParList.Dwelltimes = 2 * ParList.Dwelltimes;																													% DATASET IS NOT SPIRAL!
	else																																								% THIS IS A HACK!
		fprintf('\n\nWARNING: I DID  N O T  DOUBLE THE DWELLTIMES AS USUAL.')
		fprintf('\nIF YOUR DATASET IS A CONVENTIONAL, FULLY SAMPLED (no elliptical or acquisition weighting) DATASET,\nTHE RESULTS WILL BE WRONG!')
	end
end

if(numel(strfind(file_path, '.dat')) > 0 && ~isempty(regexpi(ParList.AssumedSequence,'Imaging_')))
    if(ParList.ThreeD_flag && ParList.nPartEnc > 1)  
        ParList.FoV_Partition = 2*ParList.FoV_Partition;
    else
        ParList.FoV_Read = 2*ParList.FoV_Read;
    end
end


if(ParList.ThreeD_flag && ParList.nPartEnc == 1)
    ParList.ThreeD_flag = false;
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


% Add field gyromagnetic ratio based on nucleus
if(~isempty(regexp(ParList.Nucleus,'1H','ONCE')))
    ParList.GyroMagnRatioOverTwoPi = 42.57747892 * 10^6;
elseif(~isempty(regexp(ParList.Nucleus,'2H','ONCE')))
    ParList.GyroMagnRatioOverTwoPi =  6.5357349913 * 10^6;
elseif(~isempty(regexp(ParList.Nucleus,'31P','ONCE')))
    ParList.GyroMagnRatioOverTwoPi =  17.25263 * 10^6;
end










% Consistency Checkings

% Check if All Vectors for Multislice Data have same size.
% It could happen that some are initialized by zeros with dimension [1 1],
% and not set because not available in the header,
% but 4 slices were measured, and some other values have dimension [1 4]. 
% This could lead to problems, e.g. the following:
% SliceNormalVector_x = 0; SliceNormalVector_y = [0.12 0.12]; SliceNormalVector_z = [0.9 0.9];

MaxSizeSliceNormalVecs = max(cat(1,numel(ParList.SliceNormalVector_x),numel(ParList.SliceNormalVector_y),numel(ParList.SliceNormalVector_z))); %#ok
%ParList.SliceNormalVector_x = repmat(ParList.SliceNormalVector_x, [1 1+MaxSizeSliceNormalVecs-numel(ParList.SliceNormalVector_x)]);

for Dim = {'x','y','z'}
    eval([ 'ParList.SliceNormalVector_' Dim{1} ' = repmat(ParList.SliceNormalVector_' Dim{1} ', [1+MaxSizeSliceNormalVecs-numel(ParList.SliceNormalVector_' Dim{1} ') 1]);' ]);
end







end

fclose(fid);


end

function ParList = CheckAssumedSequence(ParList,file_path)

    if(numel(strfind(file_path, '.dat')) > 0)
        TmpTwixHdr = mapVBVD_ReadOnlyHdr(file_path);
        ParList.SequenceDescription = TmpTwixHdr.hdr.Config.SequenceDescription;
        ParList.SequenceString = TmpTwixHdr.hdr.Config.SequenceString;
    end
    CheckFields = {ParList.tSequenceFileName};
    if(isfield(ParList,'SequenceDescription'))
        CheckFields(2) = {ParList.SequenceDescription};
    end
    if(isfield(ParList,'SequenceString'))
        CheckFields(3) = {ParList.SequenceString};
    end
    
    CheckFor = {'Eccentric','','';'CRT','Rollercoa','CONCEPT';'Spiral','','';'gre','tfl','bow_ph_map'};
    AssumedSequences = {'AntoinesEccentricOrRosette';'ViennaCRT';'BorjanSpiral';'Imaging_GRE'};
    SpatialSpectralEncoding_flags = [1;1;1;0];
    
    ParList.AssumedSequence = 'CSIOrSVS';
    ParList.SpatialSpectralEncoding_flag = 0;
    for ii = 1:size(CheckFor,1)
        for jj = 1:size(CheckFor,2)
            for kk = 1:numel(CheckFields)
                if(~isempty(CheckFor{ii,jj}) && ~isempty(CheckFields{kk}))
                    if(contains(CheckFields{kk},CheckFor{ii,jj},'IgnoreCase',true))
                        ParList.AssumedSequence = AssumedSequences{ii};
                        ParList.SpatialSpectralEncoding_flag = SpatialSpectralEncoding_flags(ii);
                        break;
                    end
                end
            end
        end
    end

end





