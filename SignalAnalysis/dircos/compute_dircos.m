function [Phase_new, Read_new] = compute_dircos(Slice_new_OrPars,Rot_angle_around_Slice_new)

% Assume that Slice_new_OrPars is ascconv part of header read in by my function "read_ascconv".
if(isstruct(Slice_new_OrPars))
    Rot_angle_around_Slice_new = Slice_new_OrPars.InPlaneRotation;
    Slice_new = [Slice_new_OrPars.SliceNormalVector_x(1) Slice_new_OrPars.SliceNormalVector_y(1) Slice_new_OrPars.SliceNormalVector_z(1)];
else
    Slice_new = Slice_new_OrPars;
end

Phase_new = [0 Slice_new(3) -Slice_new(2)] * 1/sqrt(Slice_new(2)^2+Slice_new(3)^2);
Read_new = cross(Slice_new,Phase_new);


if(Rot_angle_around_Slice_new ~= 0)
%     Rot_mat_around_Slice_new = [cos(Rot_angle_around_Slice_new) sin(Rot_angle_around_Slice_new) 0; -sin(Rot_angle_around_Slice_new) cos(Rot_angle_around_Slice_new) 0; 0 0 1];
%     Phase_new = Rot_mat_around_Slice_new * transpose(Phase_new);    
%     Read_new = Rot_mat_around_Slice_new * transpose(Read_new);
%     %Read_new = cross(Slice_new,Phase_new);
%     
%     Phase_new = transpose(Phase_new);
%     Read_new = transpose(Read_new);
    
    for i = 1:3
        Phase_new(i) = cos(Rot_angle_around_Slice_new)*Phase_new(i) - sin(Rot_angle_around_Slice_new)*Read_new(i);
    end
    
    Read_new = cross(Slice_new,Phase_new);

end

