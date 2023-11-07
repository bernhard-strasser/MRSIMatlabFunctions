function [DataStruct] = io_ReadEccentricTraj(DataStruct,TrajFile,Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations
if(~exist('Settings','var'))
    Settings = [];
end
if(~isfield(DataStruct.Par,'GradDelay_x_us'))
    if(isfield(Settings,'GradDelay_x_us'))
        DataStruct.Par.GradDelay_x_us = Settings.GradDelay_x_us;
    else
        DataStruct.Par.GradDelay_x_us = 0;
    end
end
if(~isfield(DataStruct.Par,'GradDelay_y_us'))
    if(isfield(Settings,'GradDelay_x_us'))
        DataStruct.Par.GradDelay_y_us = Settings.GradDelay_y_us;
else
        DataStruct.Par.GradDelay_y_us = 0;
    end    
end

if(~isfield(Settings,'IncludeRewinder_flag'))
    Settings.IncludeRewinder_flag = false;
end
if(~isfield(Settings,'maxR'))
    Settings.maxR = 0.5;
end
if(~isfield(Settings,'ShowTrajs'))
    Settings.ShowTrajs = false;
end



DataStruct.Par.GradDelay_y_us = 0;
DataStruct.Par.GradDelay_x_us = 0;


DataStruct.Par.RewPts = 0;


%% Run Files 
if(~exist('TrajFile','var') || isempty(TrajFile))

%     for a=1:DataStruct.Par.nAngInts
%         if (DataStruct.Par.WipMemBlock_alFree(8)==1) % Rosette
%              LoopCenter_Y = -DataStruct.Par.LoopCenter_Rad(a)*cos(DataStruct.Par.LoopCenter_Ph(a));
%              LoopCenter_X = -DataStruct.Par.LoopCenter_Rad(a)*sin(DataStruct.Par.LoopCenter_Ph(a)); 
%              LoopRadius = DataStruct.Par.LoopRadius(a);
%              LoopStartAngleToRA = DataStruct.Par.LoopStartAngletoRA(a)-pi;
% 
%              for s=1:DataStruct.Par.TrajPts(a)
%                  
%                 DataStruct.InTraj.GM{a}(1,s) = LoopCenter_X+LoopRadius*sin(2*pi*(s-3)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );      
%                 DataStruct.InTraj.GM{a}(2,s) = LoopCenter_Y+LoopRadius*cos(2*pi*(s-3)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );                  
% %             DataStruct.InTraj.GM{a}(:,s) = [LoopRadius*real(exp(1j*(-s+3)*2*pi/max(DataStruct.Par.TrajPts(a)) + 1j*LoopStartAngleToRA )) + LoopCenter_X, LoopRadius*imag(exp(1j*(-s+3)*2*pi/max(DataStruct.Par.TrajPts(a)) + 1j*LoopStartAngleToRA ))+ LoopCenter_Y];
%              
%              end
%         elseif (DataStruct.Par.WipMemBlock_alFree(8)==3) % Ellipse
%              LoopRadius = DataStruct.Par.LoopRadius(a);
%              LoopStartAngleToRA = DataStruct.Par.LoopStartAngletoRA(1)-pi/2;
%              LoopCenter_Ph = DataStruct.Par.LoopCenter_Ph(a);
%              
%              R = [cos(LoopCenter_Ph) -sin(LoopCenter_Ph); sin(LoopCenter_Ph) cos(LoopCenter_Ph)];
% %              Efakt = DataStruct.Par.WipMemBlock_adFree(14);
%              Efakt = DataStruct.Par.WipMemBlock_adFree(8);
%              
%              for s=1:DataStruct.Par.TrajPts(a)
%                 
%                 DataStruct.InTraj.GM{a}(2,s) = Efakt*LoopRadius*sin(2*pi*(s-3)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );      
%                 DataStruct.InTraj.GM{a}(1,s) = -LoopRadius-LoopRadius*cos(2*pi*(s-3)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  
%                 DataStruct.InTraj.GM{a}(:,s) = DataStruct.InTraj.GM{a}(:,s)'*R;
%                 
%              end
%          elseif (DataStruct.Par.WipMemBlock_alFree(8)==4) % EGG 
%              LoopRadius = DataStruct.Par.LoopRadius(a);
%              LoopStartAngleToRA = DataStruct.Par.LoopStartAngletoRA(1)-pi/2;
%              LoopCenter_Ph = DataStruct.Par.LoopCenter_Ph(a);
%              
%              R = [cos(LoopCenter_Ph) -sin(LoopCenter_Ph); sin(LoopCenter_Ph) cos(LoopCenter_Ph)];
%              Efakt = DataStruct.Par.WipMemBlock_adFree(8);
%              Efakt_2 = DataStruct.Par.WipMemBlock_adFree(9);
%              for s=1:DataStruct.Par.TrajPts(a)
%                 
%                 DataStruct.InTraj.GM{a}(2,s) = (Efakt_2-(sin(pi*(s-3)/max(DataStruct.Par.TrajPts(a))))*(Efakt_2-Efakt))*LoopRadius*sin(2*pi*(s-3)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );
%                 DataStruct.InTraj.GM{a}(1,s) = -LoopRadius-LoopRadius*cos(2*pi*(s-3)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  
%                 DataStruct.InTraj.GM{a}(:,s) = DataStruct.InTraj.GM{a}(:,s)'*R;
%                 
%              end        
%          else   
%                 fprintf('\n wrong trajectory type, please select different reconstruction')
%                 break
%         end
%     end

    ShiftTrajPt = 2;
    for a=1:DataStruct.Par.nAngInts
        if (DataStruct.Par.WipMemBlock_alFree(8)==1) % Rosette
             LoopCenter_Y = -DataStruct.Par.LoopCenter_Rad(a)*cos(DataStruct.Par.LoopCenter_Ph(a));
             LoopCenter_X = -DataStruct.Par.LoopCenter_Rad(a)*sin(DataStruct.Par.LoopCenter_Ph(a)); 
             LoopRadius = DataStruct.Par.LoopRadius(a);
             LoopStartAngleToRA = DataStruct.Par.LoopStartAngletoRA(a)-pi;

             for s=1:DataStruct.Par.TrajPts(a)
                 
                DataStruct.InTraj.GM{a}(1,s) = LoopCenter_X+LoopRadius*sin(2*pi*(s-ShiftTrajPt)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  % bstr: changed from -3 to -2    
                DataStruct.InTraj.GM{a}(2,s) = LoopCenter_Y+LoopRadius*cos(2*pi*(s-ShiftTrajPt)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  % bstr: changed from -3 to -2                  
%             DataStruct.InTraj.GM{a}(:,s) = [LoopRadius*real(exp(1j*(-s+3)*2*pi/max(DataStruct.Par.TrajPts(a)) + 1j*LoopStartAngleToRA )) + LoopCenter_X, LoopRadius*imag(exp(1j*(-s+3)*2*pi/max(DataStruct.Par.TrajPts(a)) + 1j*LoopStartAngleToRA ))+ LoopCenter_Y];
             
             end
        elseif (DataStruct.Par.WipMemBlock_alFree(8)==3) % Ellipse
             LoopRadius = DataStruct.Par.LoopRadius(a);
             LoopStartAngleToRA = DataStruct.Par.LoopStartAngletoRA(1)-pi/2;
             LoopCenter_Ph = DataStruct.Par.LoopCenter_Ph(a);
             
             R = [cos(LoopCenter_Ph) -sin(LoopCenter_Ph); sin(LoopCenter_Ph) cos(LoopCenter_Ph)];
%              Efakt = DataStruct.Par.WipMemBlock_adFree(14);
             Efakt = DataStruct.Par.WipMemBlock_adFree(8);
             
             for s=1:DataStruct.Par.TrajPts(a)
                
                DataStruct.InTraj.GM{a}(2,s) = Efakt*LoopRadius*sin(2*pi*(s-ShiftTrajPt)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  % bstr: changed from -3 to -2      
                DataStruct.InTraj.GM{a}(1,s) = -LoopRadius-LoopRadius*cos(2*pi*(s-ShiftTrajPt)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  % bstr: changed from -3 to -2  
                DataStruct.InTraj.GM{a}(:,s) = DataStruct.InTraj.GM{a}(:,s)'*R;
                
             end
         elseif (DataStruct.Par.WipMemBlock_alFree(8)==4) % EGG 
             LoopRadius = DataStruct.Par.LoopRadius(a);
             LoopStartAngleToRA = DataStruct.Par.LoopStartAngletoRA(1)-pi/2;
             LoopCenter_Ph = DataStruct.Par.LoopCenter_Ph(a);
%              fact_2 = 2.8;
             
             R = [cos(LoopCenter_Ph) -sin(LoopCenter_Ph); sin(LoopCenter_Ph) cos(LoopCenter_Ph)];
             Efakt = DataStruct.Par.WipMemBlock_adFree(8);
             Efakt_2 = DataStruct.Par.WipMemBlock_adFree(9);
             for s=1:DataStruct.Par.TrajPts(a)
%                 fact = (fact_2*1/Efakt_2-(sin(pi*(s)/max(DataStruct.Par.TrajPts(a))))*(fact_2*1/Efakt_2-fact_2*1/Efakt));
                CorrVal1 = DataStruct.Par.WipMemBlock_adFree(8)+1;
                CorrVal2 = DataStruct.Par.WipMemBlock_adFree(9)+1;
                fact = Efakt+1-(sin(pi*(s-ShiftTrajPt)/max(DataStruct.Par.TrajPts(a))))*(Efakt-Efakt_2);
                fact = 0;
                ShiftTrajPt2 = 0;
%                 ShiftTrajPt2 = ShiftTrajPt;
                
                DataStruct.InTraj.GM{a}(2,s) = (Efakt_2-(sin(pi*(s-ShiftTrajPt2)/max(DataStruct.Par.TrajPts(a))))*(Efakt_2-Efakt))*LoopRadius*sin(2*pi*(s-fact)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );
                DataStruct.InTraj.GM{a}(1,s) = -LoopRadius-LoopRadius*cos(2*pi*(s-ShiftTrajPt2)/max(DataStruct.Par.TrajPts(a)) + LoopStartAngleToRA );  
                DataStruct.InTraj.GM{a}(:,s) = DataStruct.InTraj.GM{a}(:,s)'*R;
                
                
             end    
             DataStruct.InTraj.GM{a} = circshift(DataStruct.InTraj.GM{a},[0 ShiftTrajPt]);  % I dont know how to shift the trajectory otherwise. Above code does not just 
                                                                                            % shift data.

         else   
                fprintf('\n wrong trajectory type, please select different reconstruction')
                break
        end
    end
    
    DataStruct.InTraj.GM = cellfun(@(x) x*1000, DataStruct.InTraj.GM,'uni',false);
    
else
if(~iscell(TrajFile) && ~endsWith(TrajFile,'.m'))
    Files = dir(TrajFile);
    Files = {Files.name};
    Files = Files(endsWith(Files,'.m'));
    TrajFile = strcat(TrajFile,'/',Files);
end
dGradientValues = cell(1);
for ii = 1:numel(TrajFile)
    run(TrajFile{ii});
    CurLoop = numel(dGradientValues);                       % The files are not ordered correctly, so need this complication
    dGradientValues2(CurLoop) = dGradientValues(CurLoop);
    clear dGradientValues
    NumberOfBrakeRunPointsCell{CurLoop} = NumberOfBrakeRunPoints(CurLoop);
    NumberOfLaunTrackPointsCell{CurLoop} = NumberOfLaunTrackPoints(CurLoop);
    NumberOfLoopPointsCell{CurLoop} = NumberOfLoopPoints(CurLoop);
end
clear Files TrajFile
dGradientValues = dGradientValues2; clear dGradientValues2;

    % SB2022: numer of Trajectories is smaller for Ellipses 
    DataStruct.Par.nAngInts=numel(dGradientValues);




%% Interpolate to ADC

% The points in our trajectory are GradRasterTime apart from each other.
% Simulate ADC_dt = GradientRaterTime/2

for ii = 1:numel(dGradientValues)
        
    
    % Calculate Gradient Moments of LaunchTracks
    GMLaunchtrack = -trapz(dGradientValues{ii}(1:NumberOfLaunTrackPointsCell{ii})*dMaxGradAmpl*10);
    
    
    
    CurTraj = dGradientValues{ii}(NumberOfLaunTrackPointsCell{ii}:NumberOfLaunTrackPointsCell{ii}+NumberOfLoopPointsCell{ii}-1);    

    % The trajectory is circular, so append the first 3 points to the end, and the last ones to the beginning for better interpolation
    % Append for now always the whole trajectory in beginning and end. When measuring few ADC points for the inner circles, I actually have
    % ADC-dt > GRAD_RASTER_TIME. Now the 6 extra gradient points are not always dividable by ADC-dt, but if I append the whole trajectory, it is.
    AppendPtsInBeginAndEnd = numel(CurTraj);
    CurTraj = cat(2,CurTraj(end-AppendPtsInBeginAndEnd+1:end),CurTraj,CurTraj(1:AppendPtsInBeginAndEnd));
    
    % Interpolate to ADC-grid
    DataStruct.Par.ADC_OverSamp = 1E4/DataStruct.Par.ADCdtPerAngInt_ns(ii);
    CurTraj = interp1(1:1:numel(CurTraj),CurTraj,1:(1/DataStruct.Par.ADC_OverSamp):numel(CurTraj),'spline');

    % Now remove the extra points again at beginning and end
    TakePtsFrom = DataStruct.Par.ADC_OverSamp*AppendPtsInBeginAndEnd+1; % Lets say we interpolate to 3 x finer grid. We added 3 points, totalling 3*3=9 additional points. Should start from point 10.
    TakePtsTo = TakePtsFrom + NumberOfLoopPointsCell{ii}*DataStruct.Par.ADC_OverSamp-1;   % Start from TakePtsFrom. Then we take all interpolated points from there, but minus 1
    CurTraj = CurTraj(TakePtsFrom:TakePtsTo);                               % (this last point is actually the first point of the next circumference of the circle)

    blaa = -cumtrapz(CurTraj*dMaxGradAmpl*10/DataStruct.Par.ADC_OverSamp);              % 5 is the ADC_dwelltime in us, GradientRaterTime = 2*ADC_dwelltime
    blaa = GMLaunchtrack + blaa;                                                   % Add Gradient Moment of LaunchTrack
    DataStruct.InTraj.GM{ii}(:,:) = [imag(blaa); real(blaa)];                % For some reason the x- and y-axes of the data seem to be swapped wrt to e.g. spirals...
    DataStruct.InTraj.GV{ii}(:,:) = [imag(CurTraj); real(CurTraj)];      
%     DataStruct.InTraj.GM(:,:,ii) = [real(blaa); imag(blaa)];                % For some reason the x- and y-axes of the data seem to be swapped wrt to e.g. spirals...
%     DataStruct.InTraj.GV(:,:,ii) = [real(CurTraj); imag(CurTraj)]; 
end
end 
%% Normalize DataStruct.InTraj

DataStruct.InTraj.GM = cellfun(@(x) x/(Settings.maxR*2), DataStruct.InTraj.GM,'uni',false);
% DataStruct.InTraj.GM = cellfun(@(x) x/10, DataStruct.InTraj.GM,'uni',false);



%% Debug: Show Trajectories
if(Settings.ShowTrajs)
    figure;
    if(isfield(DataStruct,'OutTraj'))
        scatter(squeeze(DataStruct.OutTraj.GM(1,1,:)),squeeze(DataStruct.OutTraj.GM(2,1,:)),'b');
    end
    
    hold on
    for AngIntNo = 1:DataStruct.Par.nAngInts
        scatter(squeeze(DataStruct.InTraj.GM{AngIntNo}(1,:)), squeeze(DataStruct.InTraj.GM{AngIntNo}(2,:)),'r')
        plot(squeeze(DataStruct.InTraj.GM{AngIntNo}(1,:)), squeeze(DataStruct.InTraj.GM{AngIntNo}(2,:)),'r')
    end
    hold off
end


%% Postparations

if(~isfield(DataStruct.Par,'TrajPts'))
    DataStruct.Par.TrajPts = cellfun(@(x) size(x,2) ,DataStruct.InTraj.GM);
end
DataStruct = supp_UpdateRecoSteps(DataStruct,Settings);

