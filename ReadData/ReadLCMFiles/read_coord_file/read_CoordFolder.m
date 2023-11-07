function Out = read_CoordFolder(file_path,MRSISizeOrMRSIPar)
%
% read_coord_file Read All Coordinate LCModel Files Within Folder
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


if(exist('MRSISizeOrMRSIPar','var'))
    if(isstruct(MRSISizeOrMRSIPar))
        MRSISize = [];
    else
        MRSISize = MRSISizeOrMRSIPar;
    end
end

FolderWasTarFile_flag = false;
if(endsWith(file_path,'.gz'))
    FolderWasTarFile_flag = true;
    Filepartss = fileparts(file_path);
    tmpdir = [Filepartss '/tmpppp_' num2str(randi(10^6))];
    mkdir(tmpdir);
    untar(file_path,tmpdir);
    file_path = [tmpdir '/CoordFiles'];
end



%% Read Coord-Files

CoordFiles = dir(file_path);
CoordFiles = {CoordFiles.name};
CoordFiles(1:2) = [];

XCoord = regexp(cellstr(CoordFiles),'_x\d+','match'); XCoord = [XCoord{:}]; XCoord = regexp(XCoord,'\d+','match'); XCoord = [XCoord{:}]; XCoord = cellfun(@str2double,XCoord);
YCoord = regexp(cellstr(CoordFiles),'_y\d+','match'); YCoord = [YCoord{:}]; YCoord = regexp(YCoord,'\d+','match'); YCoord = [YCoord{:}]; YCoord = cellfun(@str2double,YCoord);
ZCoord = regexp(cellstr(CoordFiles),'_z\d+','match'); ZCoord = [ZCoord{:}]; ZCoord = regexp(ZCoord,'\d+','match'); ZCoord = [ZCoord{:}]; ZCoord = cellfun(@str2double,ZCoord);

% Guess MRSISize, However, this will be usually wrong, if the mask is smaller than the original CSI size
if(~exist('MRSISize','var'))
    % Extract all x,y,z values from file-name
    MRSISize(1) = max(XCoord);
    MRSISize(2) = max(YCoord);
    MRSISize(3) = max(ZCoord);
end

CoordFiles_Full = strcat(file_path,'/',CoordFiles);
FirstCoordFile = read_coord_file(CoordFiles_Full{1});
Out.SpectrumMap = zeros([MRSISize(1) MRSISize(2) MRSISize(3) FirstCoordFile.NumberOfPoints],'single');
Out.FitMap = Out.SpectrumMap;
Out.BaselineMap = Out.SpectrumMap;

for Ind = 1:numel(CoordFiles_Full)

    xx = XCoord(Ind); yy = YCoord(Ind); zz = ZCoord(Ind);
    xx_Nii = MRSISize(1) - xx + 1;
    yy_Nii = MRSISize(2) - yy + 1;


    FileName = CoordFiles_Full{Ind};
    if(~exist(FileName,'file'))
        continue
    end

    % Read Coord File
    CurCoord = read_coord_file(FileName);
    if(numel(CurCoord.ppm_points) == 1 && CurCoord.ppm_points == 0)
        continue
    end

    % Assign Data to correct position
    Out.SpectrumMap(xx_Nii,yy_Nii,zz,:) = CurCoord.Spec;
    Out.FitMap(xx_Nii,yy_Nii,zz,:) = CurCoord.Fit;
    Out.BaselineMap(xx_Nii,yy_Nii,zz,:) = CurCoord.Baseline;
    Out.PPMExample = CurCoord.ppm_points;
                
end



%% 

if(FolderWasTarFile_flag)
    rmdir(tmpdir,'s')
    
end

