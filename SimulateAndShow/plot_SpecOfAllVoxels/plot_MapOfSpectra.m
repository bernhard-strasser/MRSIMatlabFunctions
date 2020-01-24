function [] = plot_MapOfSpectra(InStruct, Settings, Mask)
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
if (~isfield(Settings,'x_start') || isempty(Settings.x_start))
    Settings.x_start = 1;
end
if (~isfield(Settings,'x_size') || isempty(Settings.x_size))
    Settings.x_size = size(InStruct.Data,1);
end
if (~isfield(Settings,'y_start') || isempty(Settings.y_start))
    Settings.y_start = 1;
end
if (~isfield(Settings,'y_size') || isempty(Settings.y_size))
    Settings.y_size = size(InStruct.Data,2);
end
if (~isfield(Settings,'f_first') || isempty(Settings.f_first))
    Settings.f_first = 1;
end
if(~isfield(Settings,'f_end') || isempty(Settings.f_end))
    Settings.f_end = size(InStruct.Data,4);
end
if(~isfield(Settings,'plot_z') || isempty(Settings.plot_z))
    Settings.plot_z = floor(size(InStruct.Data,3)/2)+1;
end
if(~isfield(Settings,'plot_linewidth') || isempty(Settings.plot_linewidth))
    Settings.plot_linewidth = 1;
end
if(~isfield(Settings,'PlotPPM_flag') || isempty(Settings.PlotPPM_flag))
    Settings.PlotPPM_flag = true;
end
if(~isfield(Settings,'UseThisInStructMask') && ~exist('Mask','var'))
    if(isfield(InStruct,'BrainMask'))
        Settings.UseThisInStructMask = 'BrainMask';
    elseif(isfield(InStruct,'Mask'))
        Settings.UseThisInStructMask = 'Mask';        
    end
end

OpenFigs = get(groot, 'Children');
if(~isempty(OpenFigs)) 
    OpenFigs = [OpenFigs(:).Number];
else
    OpenFigs = [];
end
if(~isfield(Settings,'figno_1') || isempty(Settings.Settings.figno_1))
    Settings.figno_1 = min(setdiff(1:99,OpenFigs)); % bstr: Dont overwrite existing figures
end
if(~isfield(Settings,'figno_2') || isempty(Settings.Settings.figno_2))
    OpenFigs = [OpenFigs Settings.figno_1];
    Settings.figno_2 = min(setdiff(1:99,OpenFigs)); % bstr: Dont overwrite existing figures
end

if(isfield(Settings,'UseThisInStructMask'))
    Mask = InStruct.(Settings.UseThisInStructMask);
end
if(exist('Mask','var') && ~isempty(Mask))
    InStruct.Data = InStruct.Data .* Mask;
    [Settings.x_start,Settings.y_start] = find(Mask);
    Settings.x_end = max(Settings.x_start); Settings.y_end = max(Settings.y_start);
    Settings.x_start = min(Settings.x_start); Settings.y_start = min(Settings.y_start);
    Settings.x_size = Settings.x_end - Settings.x_start + 1; Settings.y_size = Settings.y_end - Settings.y_start + 1;
end


%%
Data_fft = fftshift(fft(InStruct.Data,[],4),4);
Data_ = abs(Data_fft(Settings.x_start:Settings.x_start+Settings.x_size-1, Settings.y_start:Settings.y_start+Settings.y_size-1,Settings.plot_z,Settings.f_first:Settings.f_end));
y_max = max(Data_(:));


pixels = size(Data_fft,2) * size(Data_fft,1);
nx = sqrt(pixels);


nnx = Settings.x_size;
nny = Settings.y_size;

widthx = 1/nnx*.9;
widthy = 1/nny*.9;
marginx = 1/nnx*0.05;
marginy = 1/nny*0.05;
originx = marginx:1/nnx:1;
originy = 1-1/nny+marginy:-1/nny:0;


xorig = Settings.x_start;
yorig = Settings.y_start;


x = xorig:xorig+nnx-1;
y = yorig:yorig+nny-1;

if(Settings.PlotPPM_flag)
    if(isfield(InStruct,'RecoPar'))
        ChemyOrPtsVec = compute_chemshift_vector(InStruct.RecoPar);
    else
        ChemyOrPtsVec = compute_chemshift_vector(InStruct.Par);        
    end
else
    ChemyOrPtsVec = 1:size(Data_fft,3);
end
PlotVec = [ChemyOrPtsVec(Settings.f_first),ChemyOrPtsVec(Settings.f_end)];

figure(Settings.figno_1);
for i=1:nnx
    for j=1:nny
%        subplot('Position', [originx(i), originy(j), widthx, widthy])
        subplot('Position', [originy(nny-j+1), originx(nnx-i+1), widthy, widthx])
        
        plot(ChemyOrPtsVec(Settings.f_first:Settings.f_end), abs((squeeze(Data_fft(x(i),y(j),Settings.plot_z, Settings.f_first:Settings.f_end)))),'k', 'LineWidth', Settings.plot_linewidth); 
        axis([min(PlotVec),max(PlotVec), 0, y_max/1]); axis off;

    end
end


ref_image = sum(abs(Data_fft(:,:,Settings.plot_z,Settings.f_first:Settings.f_end)),4);
ref_image = ref_image / max(abs(ref_image(:)));

% dB scale
% ref_image = 20*log10(sum(abs(Data_fft(:,:,Settings.f_first:Settings.f_end)),3));

figure(Settings.figno_2); imagesc(abs(ref_image)); colormap jet; colorbar; axis square; hold on;
rectangle('Position',[y(1),x(1),y(end)-y(1),x(end)-x(1)],'EdgeColor','y', 'LineStyle','--');
hold off;

