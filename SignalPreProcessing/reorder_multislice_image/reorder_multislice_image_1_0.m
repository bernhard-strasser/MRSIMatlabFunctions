function mat_out = reorder_multislice_image_1_0(mat_in,reorder_direction)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%       FUNCTION TO REORDER MULTISLICE IMAGES       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fredir = frequency encoding direction; phadir = phase encoding direction




%% 0. Preparations

reorder_size = size(mat_in, reorder_direction);
reorder_step = ceil(reorder_size/2);
circshift_number = floor(reorder_size/2);
circshift_vec = zeros([1 numel(size(mat_in))]);
circshift_vec(reorder_direction) = -circshift_number;


%% 1. Reorder Progression

reorder_progression = ones([1 reorder_size]);                                                   % deutsch(progression) = Folge
for reorder_index = 2:reorder_size
    if(mod(reorder_index,2) == 0)
        reorder_progression(reorder_index) = reorder_progression(reorder_index-1) + reorder_step;
    else
        reorder_progression(reorder_index) = reorder_progression(reorder_index-1) - (reorder_step-1);
    end
end



%% 2. Reorder String

reorder_string = '';
for mat_sizeindex = 1:numel(size(mat_in))
   if(mat_sizeindex == reorder_direction)
       reorder_string = [reorder_string '[' sprintf('%d ',reorder_progression) '],' ];
   else
       reorder_string = [reorder_string ':,'];
   end
end
reorder_string(end) = [];




%% 2. Circshift & Reorder

mat_out = circshift(mat_in,circshift_vec);

eval([ 'mat_out = mat_out(' reorder_string ');' ]); 
