function ErrorOccurred = StretchDicomFiles(outputfolder,StretchBy,varargin)
%
% template_1_0 Do nothing specific
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         A                           ...     This is output A
% -         B                           ...     This is output B
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations
ErrorOccurred = true;





if(~exist('StretchBy','var') || (exist('StretchBy','var') && (isempty(StretchBy) || (numel(StretchBy) == 1 && StretchBy == 0)) ) )
    StretchBy = [2 4 8];
end


if(numel(varargin) == 1 && exist(varargin{1},'dir'))
	% If its a dir
	DicomFileList = dir(varargin{1});
	DicomFileList = {DicomFileList.name};
	DicomFileList(cellfun('isempty',regexp(DicomFileList,'.*\.IMA','ONCE'))) = [];
	DicomFileList = strcat(varargin{1},'/',DicomFileList);
else
	% If its not a dir, but a list of files
	DicomFileList = varargin;
	
end


if(~exist('outputfolder','var') || (exist('outputfolder','var') && (isempty(outputfolder) || (numel(outputfolder) == 1 && outputfolder == 0)) ) )
    outputfolder = regexp(DicomFileList{1},'.*/','match'); outputfolder = outputfolder{1};
	if(strcmp(outputfolder(end),'/')); outputfolder(end) = []; end
end


% 0.3 Definitions
    






%% 1. Read in files, Multiply, Write Out

CurFileInd = 0;
for CurFile = DicomFileList
	CurFileInd = CurFileInd+1;
	fprintf('\nStretching file %d of %d ...',CurFileInd,numel(DicomFileList))
	
	CurFileStr = CurFile{:};
	CurFileWOPath = regexprep(CurFile,'.*/','');
	
	CurData = dicomread(CurFileStr);
	CurHead = dicominfo(CurFileStr);
	
	for StretchByInd = 1:numel(StretchBy)
		CurOutDir = [outputfolder '/StretchedBy' num2str(StretchBy(StretchByInd))];
		if(~exist(CurOutDir,'dir'))
			mkdir(CurOutDir)
		end
		CurData = StretchBy(StretchByInd) * CurData;
		CurData(CurData<0) = 0;
		CurData(CurData>4095) = 4095;
		
		NewFilePath = [CurOutDir '/' regexprep(CurFileWOPath{:},'.IMA',['_StretchedBy' num2str(StretchBy(StretchByInd)) '.IMA'])];
		dicomwrite(CurData,NewFilePath,CurHead,'CreateMode','Copy')
		
	end
	
end





%% 3. Postparations

fprintf('\n')
% fclose(fid)

ErrorOccurred = false;




