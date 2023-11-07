function [] = plot_FIDs(InStruct, Settings, Mask)
% InStruct.Data : input spectra
% Settings.x_start, Settings.y_start : upper left corner of ROI
% Settings.x_size, Settings.y_size : size of ROI to show
% Settings.f_first, Settings.f_end : indices of first and last freqs to show, optional
% Settings.figno_1, Settings.figno_2 : figure numbers, optional

% Based on Berkin Bilgic's plot_spectra function. Modified by Bernhard Strasser. 

%% Preparations

if(~exist('Settings','var'))
    Settings = struct();
end
if(~isfield(Settings,'plot_linewidth') || isempty(Settings.plot_linewidth))
    Settings.plot_linewidth = 1;
end

if(exist('Mask','var'))
    [Tmp_x,Tmp_y,Tmp_z] = find(Mask);
    Settings.PlotVoxels = cat(2,Tmp_x,Tmp_y,Tmp_z);
end

if(~isfield(Settings,'TakeRealAbsImagComplex'))
    Settings.TakeRealAbsImagComplex = @abs;
end

if(~isfield(Settings,'PlotVoxels'))
    Settings.PlotVoxels = [1 1 1];
end
if(size(Settings.PlotVoxels,2) ~= 3)
    Settings.PlotVoxels = transpose(Settings.PlotVoxels);
end
if(size(Settings.PlotVoxels,2) ~= 3)
    error('Size of Settings.PlotVoxels must be 3 in one dimension.')
end
if(~isfield(Settings,'TakeRealAbsImagComplex'))
    Settings.TakeRealAbsImagComplex = @abs;
end


%%

TimeVec = (0:(InStruct.RecoPar.vecSize-1)) * InStruct.RecoPar.Dwelltimes(1)/1E9;


% convert from Time to FID pts
if(isfield(Settings,'PlotTimeRange'))
    Settings.t_first = FindClosestIndex(TimeVec,min(Settings.PlotTimeRange)); Settings.t_first = Settings.t_first{1};
    Settings.t_end = FindClosestIndex(TimeVec,max(Settings.PlotTimeRange)); Settings.t_end = Settings.t_end{1};
else
    Settings.t_first = 1;
    Settings.t_end = size(InStruct.Data,4);
end


NoOfVox = size(Settings.PlotVoxels,1);

Data_ = [];
for ii =1:NoOfVox
    Data_ = cat(1,Data_,InStruct.Data(Settings.PlotVoxels(ii,1),Settings.PlotVoxels(ii,2),Settings.PlotVoxels(ii,3),Settings.t_first:Settings.t_end));
end
Data_ = feval(Settings.TakeRealAbsImagComplex,Data_);
y_max = max(Data_(:));
y_min = min(Data_(:));
if(y_max == y_min)
    y_max = y_min + 0.1;
end

PlotVec = [TimeVec(Settings.t_first),TimeVec(Settings.t_end)];



SubPlotSize = ceil(sqrt(NoOfVox));

figure;
for ii=1:NoOfVox
    subplot(SubPlotSize,SubPlotSize,ii)
    
    xx = Settings.PlotVoxels(ii,1); yy = Settings.PlotVoxels(ii,2); zz = Settings.PlotVoxels(ii,3);
    
    plot(TimeVec(Settings.t_first:Settings.t_end), feval(Settings.TakeRealAbsImagComplex,squeeze(InStruct.Data(xx,yy,zz, Settings.t_first:Settings.t_end))),'k', 'LineWidth', Settings.plot_linewidth); 
    axis([min(PlotVec),max(PlotVec), y_min, y_max/1]); 
    title(['x' num2str(xx) ' y' num2str(yy) ' z' num2str(zz)])

end


