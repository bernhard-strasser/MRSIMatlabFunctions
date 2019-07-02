function equation_system = compute_dircos_equation_system(dircos)

global x3_ y3_ z3_

x1_ = dircos(1);
x2_ = dircos(2);
y1_ = dircos(3);
y2_ = dircos(4);
z1_ = dircos(5);
z2_ = dircos(6);




% global z1_ z2_ z3_
% 
% x1_ = dircos(1);
% x2_ = dircos(2);
% x3_ = dircos(3);
% y1_ = dircos(4);
% y2_ = dircos(5);
% y3_ = dircos(6);





equation_system = ... % THE NON-ROTATED VECTORS ARE NORMALIZED IN THE NEW BASIS (x1_ = <x_,x>, y1_ = <y_,x>, z1_= <z_,x> --> [x1_, y1_, z1_] = (x) in the basis '
                   [x1_^2 + y1_^2 + z1_^2 - 1; ...    
                   x2_^2 + y2_^2 + z2_^2 - 1; ...
                   ...
                   ... % THE ROTATED VECTORS ARE NORMALIZED
                   x1_^2 + x2_^2 + x3_^2 - 1; ...
                   y1_^2 + y2_^2 + y3_^2 - 1; ...
                   z1_^2 + z2_^2 + z3_^2 - 1; ...
                   ...
                   ... % THE NON-ROTATED VECTORS ARE ORTHOGONAL IN THE '-BASIS ( <x,y> = <x,x_>*<y,x_> + <x,y_>*<y,y_> + <x,z_>*<y,z_> = 0 etc.
                   x1_*x2_ + y1_*y2_ + z1_*z2_; ...
                   x1_*x3_ + y1_*y3_ + z1_*z3_; ...
                   x2_*x3_ + y2_*y3_ + z2_*z3_; ...
                   ...
                   ... % THE ROTATED VECTORS ARE ORTHOGONAL
                   x1_*y1_ + x2_*y2_ + x3_*y3_; ...
                   x1_*z1_ + x2_*z2_ + x3_*z3_; ...
                   y1_*z1_ + y2_*z2_ + y3_*z3_; ...
                   ...
                   ... % THE z'-VECTOR IS THE CROSS PRODUCT OF x' and y'
                   x2_*y3_ - x3_*y2_ - z1_; ...
                   x3_*y1_ - x1_*y3_ - z2_; ...
                   x1_*y2_ - x2_*y1_ - z3_...
                   ];