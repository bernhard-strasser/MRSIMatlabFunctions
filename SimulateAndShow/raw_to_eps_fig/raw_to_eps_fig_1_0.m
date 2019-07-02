function raw_to_eps_fig_1_0(in_dir,raw_precision,ROW,COL)
%
% 

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    

% 0.3 Preparations



if(~exist('raw_precision','var'))
    raw_precision = 'float32';
end

if(~exist('ROW','var'))
    ROW = 64;
end

if(~exist('COL','var'))
    COL = 64;
end


mkdir(sprintf('%s/EPS',in_dir))
mkdir(sprintf('%s/FIG',in_dir))
mkdir(sprintf('%s/JPG',in_dir))


pause on




%% 1. Find all .raw files


file_list = dir(sprintf('%s/*.raw', in_dir));
file_list = cellstr(char(file_list.name));



%% 2. Read in NAA map and resampled map


NAA_map_logical = ~cellfun(@isempty,strfind(file_list,'NAA_amp_map.raw'));
NAA_map_res_logical = ~cellfun(@isempty,strfind(file_list,'NAA_amp_map_res.raw'));


NAA_map = read_RawFiles_1_0([in_dir '/' file_list{NAA_map_logical}],raw_precision,ROW,COL,1);
NAA_map_res = read_RawFiles_1_0([in_dir '/' file_list{NAA_map_res_logical}],raw_precision,128,128,1);




%% 2. Loop over all .raw files

for file_no = 1:numel(file_list)
    
    file = file_list{file_no};
    filename = strrep(file,'.raw','')
    
    if(strfind(file,'_res') > 0)
        ROW = 128; COL = 128;
        ratio_to_map = NAA_map_res;
    else
        ratio_to_map = NAA_map;
    end
    
    map = read_RawFiles_1_0([in_dir '/' file],raw_precision,ROW,COL,1);

    
    
    
    if(numel(strfind(file,'amp')) > 0 && ne(file_no,find(NAA_map_logical)) && ne(file_no, find(NAA_map_logical)))
        map_ratio = map ./ ratio_to_map;
        
        
        map_rat_fig = figure('visible','off');
        imagesc(map_ratio)
        colorbar('FontSize', 26)        
        saveas(map_rat_fig,sprintf('%s/EPS/%s_RatToNAA',in_dir,filename),'epsc2')
        saveas(map_rat_fig,sprintf('%s/FIG/%s_RatToNAA',in_dir,filename),'fig')
        saveas(map_rat_fig,sprintf('%s/JPG/%s_RatToNAA',in_dir,filename),'jpg')
        close(map_rat_fig)
        
    end
    
    map_reshaped = reshape(map,[ROW*COL 1]);
    mad_of_map = mad(map_reshaped(ne(map_reshaped,0)),1);
    median_of_map = median(map_reshaped(ne(map_reshaped,0)));
    map(abs(map) > median_of_map + 4*mad_of_map) = NaN;
    
    map_fig = figure('visible','off');
    imagesc(map)
    colorbar('FontSize', 26)
    saveas(map_fig,sprintf('%s/EPS/%s',in_dir,filename),'epsc2')
    saveas(map_fig,sprintf('%s/FIG/%s',in_dir,filename),'fig')
    saveas(map_fig,sprintf('%s/JPG/%s',in_dir,filename),'jpg')
    close(map_fig)
    
    
    ROW = 64; COL = 64;
end



%% 7. THE END

pause off







