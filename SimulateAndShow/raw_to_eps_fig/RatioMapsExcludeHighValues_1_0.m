function RatioMapsExcludeHighValues_1_0(in_dir,out_dir,raw_precision,ROW,COL)
%
% 

%% 0. Declarations, Preparations, Definitions

% 0.1 Declarations


% 0.2 Definitions
    

% 0.3 Preparations

run ./dir_file.m


if(~exist('raw_precision','var'))
    raw_precision = 'float';
end

if(~exist('ROW','var'))
    ROW = 64;
end

if(~exist('COL','var'))
    COL = 64;
end




pause on




%% 1. Find all .raw files


file_list = dir(sprintf('%s/*.raw', in_dir));
file_list = cellstr(char(file_list.name));



%% 2. Read in NAA map and resampled map


NAA_map_logical = ~cellfun(@isempty,strfind(file_list,'NAA_amp_map.raw'));
NAA_map_res_logical = ~cellfun(@isempty,strfind(file_list,'NAA_amp_map_res.raw'));


NAA_map = read_RawFiles_1_0([in_dir '/' file_list{NAA_map_logical}],raw_precision,ROW,COL,1);
NAA_map_res = read_RawFiles_1_0([in_dir '/' file_list{NAA_map_res_logical}],raw_precision,128,128,1);




%% 3. Remove too high NAA values


NAA_reshaped = reshape(NAA_map,[ROW*COL 1]);
std_NAA = std(NAA_reshaped(ne(NAA_reshaped,0)),1);
mean_NAA = mean(NAA_reshaped(ne(NAA_reshaped,0)));
NAA_map(abs(NAA_map) > mean_NAA + 4*std_NAA) = NaN;

NAAres_reshaped = reshape(NAAres_map,[ROW*COL 1]);
std_NAAres = std(NAAres_reshaped(ne(NAAres_reshaped,0)),1);
mean_NAAres = mean(NAAres_reshaped(ne(NAAres_reshaped,0)));
NAAres_map(abs(NAAres_map) > mean_NAAres + 4*std_NAAres) = NaN;



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
    
    
    % Kill all outliers
    if(numel(strfind(file,'amp')) > 0)
        map_reshaped = reshape(map_map,[ROW*COL 1]);
        std_map = std(map_reshaped(ne(map_reshaped,0)),1);
        mean_map = mean(map_reshaped(ne(map_reshaped,0)));
        map(abs(map) > mean_map + 4*std_map) = NaN;
    
        
        % Create Ratio Map
        if(ne(file_no,find(NAA_map_logical)) && ne(file_no, find(NAA_map_logical)))
            map_ratio = map ./ ratio_to_map;
        end
        
        
        % CRLB Maps --> Restrict to values 0 < CRLB < 60
    else
        map(map > 60) = NaN;       
    end
    
    
    
    % Write files to raw_file
    map_fid = fopen(sprintf('%s/%s_clipped.raw',out_dir,filename),'w');
    fwrite(map_fid,map,'float');
    fclose(map_fid);
    
    if(exist('map_ratio','var')
        map_ratio_fid = fopen(sprintf('%s/%s_RatToNAA.raw',out_dir,filename),'w');
        fwrite(map_ratio_fid,map,'float');
        fclose(map_ratio_fid);       
    end
    
    
    % Set variables back
    clear map_ratio
    ROW = 64; COL = 64;
end



%% 7. THE END

pause off







