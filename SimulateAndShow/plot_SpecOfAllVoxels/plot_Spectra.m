function [] = plot_Spectra(InStruct, Settings, Mask)
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

if(~isstruct(Settings) && isnumeric(Settings))
    Tmp.PlotVoxels = Settings;
    Settings = Tmp; clear Tmp
end


if(~isfield(Settings,'plot_linewidth') || isempty(Settings.plot_linewidth))
    Settings.plot_linewidth = 1;
end
if(~isfield(Settings,'PlotPPM_flag') || isempty(Settings.PlotPPM_flag))
    Settings.PlotPPM_flag = true;
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
Data_fft = fftshift(fft(InStruct.Data,[],4),4);


if(isfield(InStruct,'RecoPar'))
    ChemyOrPtsVec2 = compute_chemshift_vector(InStruct.RecoPar);
else
    ChemyOrPtsVec2 = compute_chemshift_vector(InStruct.Par);        
end
if(Settings.PlotPPM_flag)
    ChemyOrPtsVec = ChemyOrPtsVec2;
else
    ChemyOrPtsVec = 1:size(Data_fft,4);
end

% convert from PPM to freq pts
if(isfield(Settings,'PlotPPMRange'))
    Settings.f_first = FindClosestIndex(ChemyOrPtsVec,max(Settings.PlotPPMRange)); Settings.f_first = Settings.f_first{1};
    Settings.f_end = FindClosestIndex(ChemyOrPtsVec,min(Settings.PlotPPMRange)); Settings.f_end = Settings.f_end{1};
else
    Settings.f_first = 1;
    Settings.f_end = size(Data_fft,4);
end


NoOfVox = size(Settings.PlotVoxels,1);

Data_ = [];
for ii =1:NoOfVox
    Data_ = cat(1,Data_,Data_fft(Settings.PlotVoxels(ii,1),Settings.PlotVoxels(ii,2),Settings.PlotVoxels(ii,3),Settings.f_first:Settings.f_end));
end
Data_ = feval(Settings.TakeRealAbsImagComplex,Data_);
y_max = max(Data_(:));
y_min = min(Data_(:));
if(y_max == y_min)
    y_max = y_min + 0.1;
end

PlotVec = [ChemyOrPtsVec(Settings.f_first),ChemyOrPtsVec(Settings.f_end)];



SubPlotSize = ceil(sqrt(NoOfVox));

figure;
for ii=1:NoOfVox
    subplot(SubPlotSize,SubPlotSize,ii)
    
    xx = Settings.PlotVoxels(ii,1); yy = Settings.PlotVoxels(ii,2); zz = Settings.PlotVoxels(ii,3);
    
    plot(ChemyOrPtsVec(Settings.f_first:Settings.f_end), feval(Settings.TakeRealAbsImagComplex,squeeze(Data_fft(xx,yy,zz, Settings.f_first:Settings.f_end))),'k', 'LineWidth', Settings.plot_linewidth); 
    axis([min(PlotVec),max(PlotVec), y_min, y_max/1]); 
    title(['x' num2str(xx) ' y' num2str(yy) ' z' num2str(zz)])

end


