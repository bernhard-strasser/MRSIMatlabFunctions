function [LCMOutput] = io_ReadTableFile(table_file)
% ReadMetConc reads the metabolite concentrations of a CSV-file created by LCModel
% Usage: 
% Input: Path of CSV-file; 
% LCMOutput: Metabolite Concentrations (array), Metabolite Names (array)
%
%
%
% ##################################################
% ########### Function to Hell and Back ############
% ##################################################
%            Written by Bernhard Strasser




%% 0. Preps


LCMOutput.Concentrations = [];
LCMOutput.CRLBs = [];				
LCMOutput.ExtraInfo = [];
LCMOutput.MetaboliteNames = {''};
LCMOutput.ExtraNames = {''};



%% Read Table File


fid_vox = fopen(table_file,'r');
if(fid_vox == -1)

    warning('In io_ReadTableFile: File %s not found', table_file)

else
    
    ConcFound = false;
    MiscFound = false;
    DiagFound = false;
    CurMetabo = 0;
        
    tline = fgetl(fid_vox);
    while ischar(tline)    
    
        tline = fgetl(fid_vox);
        if(isempty(tline))
            continue;
        end
        if(contains(tline, '$$CONC'))
            ConcFound = true;
        end
        if(contains(tline, '$$MISC'))
            MiscFound = true;
        end
        if(contains(tline, '$$DIAG'))
            DiagFound = true;
        end
                
        
        
        
        if(ConcFound && ~MiscFound && ~DiagFound)

            met_line = regexp(tline,'\s+','split');
            NoOfMetabos = str2double(met_line{2})-1;
            if(contains(tline, '$$CONC'))
                LCMOutput.CRLBs = zeros([1 NoOfMetabos]);
                LCMOutput.Concentrations = zeros([1 NoOfMetabos]);
                LCMOutput.MetaboliteNames = cell([1 NoOfMetabos]);
                continue;                
            end
            if(contains(tline, 'Conc.'))
                continue;
            end
            CurMetabo = CurMetabo + 1;
            met_line = regexp(tline,'\s+','split');
            LCMOutput.Concentrations(CurMetabo) = str2double(met_line{2});
            LCMOutput.CRLBs(CurMetabo)  = str2double(regexprep(met_line{3},'%',''));
            LCMOutput.MetaboliteNames(CurMetabo) = met_line(end-1);
            
        elseif(MiscFound && ~DiagFound)
        

            % read data for extra info line by line
            tline = fgetl(fid_vox);
            misc_line = regexp(tline,'\s+','split');

            LCMOutput.ExtraInfo(1) = str2double(misc_line{4});    %FWHM
            LCMOutput.ExtraInfo(2) = str2double(misc_line{8});    %SNR

            tline = fgetl(fid_vox);
            misc_line = regexp(tline,'=','split');
            misc_line = regexp(misc_line{2},'ppm','split');

            LCMOutput.ExtraInfo(3) = str2double(misc_line{1});    %shift

            tline = fgetl(fid_vox);
            misc_line = regexp(tline,':','split');
            misc_line = regexp(misc_line{2},'\s+','split');

            if isnan(str2double(misc_line{1}))
                LCMOutput.ExtraInfo(4) = str2double(misc_line{2});    %zero-order phase
                LCMOutput.ExtraInfo(5) = str2double(misc_line{4});    %first-order phase
            else
                LCMOutput.ExtraInfo(4) = str2double(misc_line{1});    %zero-order phase
                LCMOutput.ExtraInfo(5) = str2double(misc_line{3});    %first-order phase
            end
            
            % We are finished after that
            break;
            
        elseif(DiagFound)
            break;
        end    
    
    
    end
    fclose(fid_vox);

    LCMOutput.ExtraNames = {'FWHM','SNR','shift','0_pha','1_pha'};
    
end
    
    
    





