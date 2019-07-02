function [Parameters,dGradientValues, dGradMom, dGradSlew] = ShowSequenceArbGrads(PathToTxtFile,PlotFlag,PauseDuration,CalcCircle_flag)
%
% ShowSequenceArbGrads Do nothing specific
%
% This function was written by Bernhard Strasser, Jannuary 2015.
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% [A,B] = read_csi_dat_1_10(inputvar1,inputvar2)
%
% Input: 
% -         PathToTxtFile                   ...    Path to text file containing the arbitrary gradient values
%
% Output:
% -         A                           ...     This is output A
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



    %% 0. Preparations, Definitions


    % 0.1 Preparations

    if(~exist('PathToTxtFile','var'))
        fprintf('\nNo input file. Abort.')
        return
    end
    if(~exist('PlotFlag','var'))
        PlotFlag = false;
        PauseDuration = -5;
    end
    if(~exist('PauseDuration','var'))
        PauseDuration = -5;        
    end
    if(~exist('CalcCircle_flag','var'))
        CalcCircle_flag = false;        
    end
    
    % If input is directory, gather all trajectory files
    if(exist(PathToTxtFile,'dir'))
        bla = dir(PathToTxtFile);
        bla3 = {bla.name};
        bla3 = strcat(PathToTxtFile, '/', bla3);
        bla3(cellfun(@exist,bla3) ~= 2) = [];
        PathToTxtFile = bla3; clear bla*
                
        text = regexp(PathToTxtFile,'/[^/][0-9]*[^/]*\.','match');
        text2 = cellfun(@cell2mat,text,'UniformOutput',false);
        text3 = regexp(text2,'[0-9]*','match');
        text4 = cellfun(@cell2mat,text3,'UniformOutput',false);
        text5 = cellfun(@str2double,text4);
        [text6, sorty] = sort(text5);
        PathToTxtFile = PathToTxtFile(sorty); clear text*

    else
        PathToTxtFile = {PathToTxtFile};
    end

    % 0.3 Definitions






    %% 1. Read File
    Index = 0;
    dGradMom = cell([1 numel(PathToTxtFile)]); dGradientValues = dGradMom; dGradSlew = dGradMom;
    for CurFile = PathToTxtFile
        Index = Index+1;
        
        %Inside should be variables: dGradientValues, dMaxGradAmpl, NumberOfBrakeRunPoints, NumberOfLaunTrackPoints, NumberOfLoopPoints
        run(CurFile{1});

        Parameters{Index}.dMaxGradAmpl = dMaxGradAmpl;
        Parameters{Index}.NumberOfBrakeRunPoints = NumberOfBrakeRunPoints;
        Parameters{Index}.NumberOfLaunTrackPoints = NumberOfLaunTrackPoints;
        Parameters{Index}.NumberOfLoopPoints = NumberOfLoopPoints;
        
%         % What is that doing?
%         dGradientValues{Index}(isnan(dGradientValues{Index})) = [];
%         if(mod(numel(dGradientValues{Index}),3) == 1)
%             dGradientValues{Index}(end) = [];
%             size_Data(1) = size_Data(1)-1;
%         end
%         dGradientValues{Index} = reshape(dGradientValues{Index},[size_Data(1) size_Data(2)-2])*dMaxGradAmpl;

        dGradientValues{Index} = dGradientValues{Index}*dMaxGradAmpl;

        
        %% 2. Calculate Positions

        if(nargout > 2)
%             % Need to integrate, i.e. s(t_n) = sum(i<n, v(t_i) * Delta_t) = Delta_t * sum(i<n,v(t_i))            % really i<n or i<=n ?
%             GradPos_x = zeros([1 size_Data(2)+1]);
%             GradPos_y = zeros([1 size_Data(2)+1]);
%             for i=2:size_Data(2)+1
%                 GradPos_x(i) = GradPos_x(i-1) + (dGradientValues{Index}(1,i-1));     % mean: Trapezformel
%             end
%             GradPos_y = imag(GradPos_x);
%             GradPos_x = real(GradPos_x);
        
        
            GradPos_x = cumsum(dGradientValues{Index});
            GradPos_y = imag(GradPos_x);
            GradPos_x = real(GradPos_x);

            dGradMom{Index} = cat(1,GradPos_x,GradPos_y) * 10;     % in units of mT*us/m
            clear GradPos_x GradPos_y
        end

        
        %% 3. Calculate SlewRates

        if(nargout > 3)
            % Numeric differentiation: a(t_n) = (a(t_n) - a(t_n-1)) / Delta_t

            dGradSlew{Index} = -diff(dGradientValues{Index});        
            dGradSlew{Index} = dGradSlew{Index} / 10;                                              % in units of mT/m per us
            dGradSlew{Index} = cat(1,real(dGradSlew{Index}),imag(dGradSlew{Index}));
        end
        
        
        %% Reshape GradValues
        
        dGradientValues{Index} = cat(1,real(dGradientValues{Index}),imag(dGradientValues{Index}));
        
        

        %% Calculate center of circle:
        if(CalcCircle_flag == true)
            GradPosFindCenter = GradPos(round(size(GradPos,1)/3):4:end,:);
            mr = (GradPosFindCenter(2,2)-GradPosFindCenter(1,2))/(GradPosFindCenter(2,1)-GradPosFindCenter(1,1));
            mt = (GradPosFindCenter(3,2)-GradPosFindCenter(2,2))/(GradPosFindCenter(3,1)-GradPosFindCenter(2,1));
            CircleCenter(1) = (mr*mt*(GradPosFindCenter(3,2)-GradPosFindCenter(1,2)) + mr*(GradPosFindCenter(2,1)+GradPosFindCenter(3,1)) - mt*(GradPosFindCenter(1,1)+GradPosFindCenter(2,1)))/(2*(mr-mt));
            CircleCenter(2) = (GradPosFindCenter(1,2)+GradPosFindCenter(2,2))/2 - 1/mr*(CircleCenter(1)-(GradPosFindCenter(1,1)+GradPosFindCenter(2,1))/2);

            Radius =  sqrt(sum((GradPos(end,:) - CircleCenter).^2));

            GradValuesFindCenter = dGradientValues{Index}(round(size(dGradientValues{Index},1)/3):end,1:2);
            CircleCenter_GradValues(1) = min(GradValuesFindCenter(:,1)) + 0.5*(max(GradValuesFindCenter(:,1)) - min(GradValuesFindCenter(:,1)));
            CircleCenter_GradValues(2) = min(GradValuesFindCenter(:,2)) + 0.5*(max(GradValuesFindCenter(:,2)) - min(GradValuesFindCenter(:,2)));
            fprintf('\nCenter of Circle %d: [%f, %f], Radius = %f',Index,CircleCenter(1),CircleCenter(2), Radius)
        end

        
        
        
        %% Plot
        if(PlotFlag)
            close all
            Data.GradPos = GradPos; Data.GradValues = GradValues; Data.GradSlew = GradSlew; Data.CircleCenter = CircleCenter; Data.CircleCenter_GradValues = CircleCenter_GradValues;
            PlotSequenceArbGrads(Data,PauseDuration)
        end
    end



end


function PlotSequenceArbGrads(Data,PauseDuration)



    %% 3. Plot gradient values
    
    
    AutoPause = true;
    if(PauseDuration == 0)
        AutoPause = false;
    end
    
    
    % Time Plot

    % Grad Values
    GradFig = figure;
    hold on
    movegui(GradFig,'northwest')

    % Set Axis
    MaxValues = max(abs(Data.dGradientValues));
    Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
    axis(Plot_Lims)

    % Labels on axes
    xlabel('GradStrength x [mT/m]')
    ylabel('GradStrength y [mT/m]')
    scatter(Data.CircleCenter_GradValues(1),Data.CircleCenter_GradValues(2),'r')


    % Grad Pos
    GradPosFig = figure;
    hold on
    movegui(GradPosFig,'northeast')

    % Set Axis
    MaxValues = max(abs(Data.GradPos));
    Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
    axis(Plot_Lims)

    % Labels on axes
    xlabel('GradPos x [mT*us/m]')
    ylabel('GradPos y [mT*us/m]')
    scatter(Data.CircleCenter(1),Data.CircleCenter(2),'r')



    % Grad SlewRate
    GradSlewFig = figure;
    hold on
    movegui(GradSlewFig,'south')

    % Set Axis
    MaxValues = max(abs(Data.GradSlew));
    Plot_Lims = 1.1 * [-max(MaxValues) max(MaxValues) -max(MaxValues) max(MaxValues)];
    axis(Plot_Lims)

    % Labels on axes
    xlabel('GradSlew x [mT/m per us]')
    ylabel('GradSlew y [mT/m per us]')



    
    % Plot in a Loop
    pause on
    for i = 1:size(Data.GradPos,1)
        if(i>1)
            figure(GradFig)
            scatter(Data.dGradientValues(i-1,1),Data.dGradientValues(i-1,2),'filled','b')
        end

        figure(GradPosFig)
        scatter(Data.GradPos(i,1),Data.GradPos(i,2),'filled','b')

        %fprintf('\nRadius: %f',norm(Data.GradPos(i,:)))

        if(i>2)
            figure(GradSlewFig)
            scatter(Data.GradSlew(i-2,1),Data.GradSlew(i-2,2),'filled','b')  
            sqrt(Data.GradSlew(i-2,:) * transpose(Data.GradSlew(i-2,:)));
        end
        if(AutoPause)
            pause(PauseDuration)
        else
            pause
        end
    end
    pause
    hold off
    pause off

end

%% 3. Postparations

% fclose(fid)






