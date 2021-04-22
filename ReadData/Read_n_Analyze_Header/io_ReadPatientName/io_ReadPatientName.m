function [PatientName] = io_ReadPatientName(File)
%
% io_ReadPatientName Read patient name from Siemens twix or dicom file
%
% This function was written by Bernhard Strasser, September 2020.
%
%
% The function can read in the patient name from Siemens twix or dicom files
%
%
% [PatientName] = io_ReadPatientName(File)
%
% Input: 
% -         File                ...  The Siemens twix of dicom file.
%
% Output:
% -         PatientName         ...  The name of the patient
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ???

% Further remarks:



%% 0. Preparations

% Get first IMA or dat file if folder is provided
if(isfolder(File))
    Subfiles = dir(File); 
    Subfiles = {Subfiles.name};
	DicomFilePresent = endsWith(Subfiles,'.IMA');
    DatFilePresent = endsWith(Subfiles,'.dat');
    if(any(DicomFilePresent))
        File = strcat(File,filesep,Subfiles{find(DicomFilePresent,1)});
    elseif(any(DatFilePresent))
        File = strcat(File,filesep,Subfiles{find(DicomFilePresent,1)});
    else
       PatientName = 'UnknownName'; 
       return; 
    end
end




%% Read Patient Name


% Read Patientname
if(endsWith(File, '.IMA'))
    [bla, PatientName] = unix(['dcmdump +P "0010,0010" ' File] );
    PatientName = regexp(PatientName,'\[(?!^).*\^?(?!^).*\]','match');
    if(bla > 0)		% If dcmdump doesnt exist, search for the Patient field manually. One day this probably should be directly implemented in Matlab, without using bash via unix command
        [bla, Border1] = unix(['grep --color=''never'' -o -b -u -P -a ''\x10\x00\x10\x00'' ' File]);			% Search for byte-offset of field 10001000, which is Patient name
        Border1 = regexp(Border1,'\d+:','match'); Border1 = str2double(Border1{1}(1:end-1));				% Get the byte-offset
        [bla, Border2] = unix(['grep --color=''never'' -o -b -u -P -a ''\x10\x00\x20\x00'' ' File]);			% Search for byte-offset of field 10002000, which is the following field
        Border2 = regexp(Border2,'\d+:','match'); Border2 = str2double(Border2{1}(1:end-1));
        [bla, PatientName] = unix(['head -c ' num2str(Border2) ' ' File ' | tail -c ' num2str(Border2-Border1-8)]);	% Cut everything out between those fields
        PatientName = cellstr(PatientName);
    end

    if(isempty(PatientName))
        PatientName = cellstr('NoName');		% E.g. if dcmdump does not exist.
    end
    PatientName = regexprep(PatientName{:},{'\[','\]','\^'},{'','','_'});
else
    [bla, PatientName] = unix(['grep -a -m 1 "tPatientName" ' File]);		% Example Result: <ParamString."tPatientName"> { "Strasser^Bernhard" }
    PatientName = regexp(PatientName,'{ "(?!^).*\^?(?!^).*"  }','match');                       % (?!^).*: Any character 0 or more times, but no caret, \^?: ^ 0 or 1 times. --> { "Strasser^Bernhard" }
    PatientName = regexp(PatientName{:},'".*"','match');                                        % --> "Strasser^Bernhard"
    if(isempty(PatientName))
        PatientName = cellstr('NoName');		% E.g if dcmdump does not exist.
    end
    PatientName = regexprep(PatientName{:},{'"','\^'},{'','_'});                                % Replace " with nothing and caret with underscore. --> Strasser_Bernhard  
end

PatientName = regexprep(PatientName,' ','_');

if(isempty(PatientName) || ~isempty(regexp(PatientName,'xxxxxxx','ONCE')))     % VD/VE twix files have this problem that he patient name is not correctly set?
    PatientName = 'UnknownName';		% E.g. if dcmdump does not exist.
end




    

%% Postparations


