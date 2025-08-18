function [LCMOutput] = io_ReadAllTableFiles(folder,MRSISize)
% ReadMetConc reads the metabolite concentrations of a CSV-file created by LCModel
% Usage: 
% Input: Path of CSV-file; 
% Output: Metabolite Concentrations (array), Metabolite Names (array)
%
%
%
% ##################################################
% ########### Function to Hell and Back ############
% ##################################################
%            Written by Bernhard Strasser




%% 0. Preps


if(isstruct(MRSISize))
    if(isfield(MRSISize,'Data'))
        MRSISize = size(MRSISize.Data);
        MRSISize = MRSISize(1:3);
    end
else
    MRSISize = MRSISize(1:3);    
end


LCMOutput = struct;



%% Find All Table Files


FileList2 = dir([folder '/*.t*']);       % Quite lazy checking...
FileList = {FileList2.name};
clear FileList2


%% Read First Meaningful Table File

for ii = 1:numel(FileList)
    Tmp = io_ReadTableFile([folder '/' FileList{ii}]);
    if(numel(Tmp.Concentrations) > 0)
        
        NoOfMetabos = numel(Tmp.Concentrations);
        LCMOutput.Concentrations = zeros([MRSISize NoOfMetabos]);
        LCMOutput.CRLBs = zeros([MRSISize NoOfMetabos]);
        LCMOutput.ExtraInfo = zeros([MRSISize numel(Tmp.ExtraNames)]);
        LCMOutput.MetaboliteNames = Tmp.MetaboliteNames;
        LCMOutput.ExtraNames = Tmp.ExtraNames;
        break
        
    end
    if(ii == numel(FileList) && ~isfield(LCMOutput,'Concentrations'))
        warning('in io_ReadAllTableFiles: Could not find any meaningful LCModel table file');
    end
    
end


%% Read All Table Files


for ii = 1:numel(FileList)
    Tmp = io_ReadTableFile([folder '/' FileList{ii}]);
    
    if(numel(Tmp.Concentrations) > 0)
        
        % find x y z position of current table file
        xx = regexp(FileList{ii},'_x[0-9]+_','match'); xx = regexp(xx{1},'[0-9]+','match'); xx = str2double(xx{1});
        yy = regexp(FileList{ii},'_y[0-9]+_','match'); yy = regexp(yy{1},'[0-9]+','match'); yy = str2double(yy{1});
        zz = regexp(FileList{ii},'_z[0-9]+\.','match'); zz = regexp(zz{1},'[0-9]+','match'); zz = str2double(zz{1});
        
        LCMOutput.Concentrations(xx,yy,zz,:) = Tmp.Concentrations;
        LCMOutput.CRLBs(xx,yy,zz,:) = Tmp.CRLBs;
        LCMOutput.ExtraInfo(xx,yy,zz,:) = Tmp.ExtraInfo;
       
    end
    
end
    
    





