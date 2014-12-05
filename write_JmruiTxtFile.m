function Error = write_JmruiTxtFile(file,DataAndHeader)
%
% write_csi_dat Write jmrui text file.
%
% This function was written by Bernhard Strasser, Dec 2014.
%
%
% The function writes a jmrui txt file.
%
%
% Error = read_JmruiTxtFile(file,DataAndHeader)
%
% Input: 
% -         file                    ...     jmrui txt file.
% -         DataAndHeader           ...     Data and Header Info that should be written to file.
%
% Output:
% -         Error                  ...     false if no error, true if error occurred
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: 

% Further remarks: 





%% 0. Preparations

Error = true;

if(exist(file,'file'))
	fprintf('\nError: File\n''%s''\n exists.\n',file)
	return;
end



fid = fopen(file,'w+');

filename = regexp(file,'/','split');
filename = filename{end};


%% 1. Write header information


fprintf(fid,'jMRUI Data Textfile\n\n');
fprintf(fid,'Filename: %s\n\n', filename);
fprintf(fid,'PointsInDataset: %s\n', num2str(DataAndHeader.Header.total_points));
fprintf(fid,'DatasetsInFile: 1\n');
fprintf(fid,'SamplingInterval: %9.4E\n', DataAndHeader.Header.dwelltime*1000);
fprintf(fid,'ZeroOrderPhase: 0E0\n');
fprintf(fid,'BeginTime: 0E0\n');
fprintf(fid,'TransmitterFrequency: %8.4E\n', DataAndHeader.Header.Frequency*1000000);
fprintf(fid,'MagneticField: 0E0\n');
fprintf(fid,'TypeOfNucleus: 0E0\n');
fprintf(fid,'NameOfPatient: \n');
fprintf(fid,'DateOfExperiment: \n');
fprintf(fid,'Spectrometer: \n');
fprintf(fid,'AdditionalInfo: \n\n\n');

fprintf(fid,'Signal and FFT\n');
fprintf(fid,'sig(real)\tsig(imag)\tfft(real)\tfft(imag)\n');
fprintf(fid,'Signal 1 out of 1 in file\n');




%% 2. Write Data


        
for i = 1:DataAndHeader.Header.total_points

	fprintf(fid,'%10.4E\t%10.4E\t%10.4E\t%10.4E\n', DataAndHeader.Data(i,1),DataAndHeader.Data(i,2),DataAndHeader.Data(i,3),DataAndHeader.Data(i,4));
	
end
	
	

%% 7. Postparations

fclose(fid);

Error = false;
