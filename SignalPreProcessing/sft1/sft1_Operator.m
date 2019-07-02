function sft1_Oper = sft1_Operator(InTraj,OutTraj,Ift_flag,InputSize)
%
% sft2 Perform 2-dimensional slow Fourier transform
%
% This function was written by Bernhard Strasser, April 2018.
%
%
% The function performs a slow Fourier transform in two spatial dimensions by using the definition of the Fourier Transform
% FT(f)(a_j) = sum(k=-N/2;N/2-1) f(x_k)exp(-2pi i x_k a_j)
% (For even N. For odd, change sum slightly)
%
%
% sft2Operator = sft2(A,Ifft_flag)
%
% Input: 
% -         Ift_flag                    ...     Flag of 0 or 1 determining if an inverse Fourier transform or normal should be performed.
% -         InputSize                           ...     Input data.
%
% Output:
% -         sft2Operator                      ...     Fourier transformation matrix which can be applied to the data by Out = sft2Operator*In
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ?







%% 0. Preparations


% Define equidistant raster so that user doesn't have to provide it
if(exist('InputSize','var'))
    % Input Sizes
    N1 = InputSize(1);
    N1_floor = floor(N1/2);

    if(~Ift_flag)          
        InScale = N1;                 % If doing an FT, assume InTraj is k-Space trajectory --> Scale k-space traj from [-0.5,0.5)
        OutScale = 1;                  % And OutTraj is i-Space trajectory --> Scale i-space traj from [-N/2, N/2)
        InTrajShift = 0;
        OutTrajShift = N1_floor;        
    else
        InScale = 1;
        OutScale = N1;
        InTrajShift = N1_floor;
        OutTrajShift = 0;
    end

    if(Ift_flag)
        
    else
        
    end
    
    if(~exist('InTraj','var') || isempty(InTraj))
        InTraj = ((1:N1) - InTrajShift - 1)/InScale;
        InTraj(isnan(InTraj)) = 0;
        clear InTraj_x InTraj_y
    end
    if(~exist('OutTraj','var') || isempty(OutTraj))
        OutTraj = ((1:N1) - OutTrajShift - 1)/OutScale;
        OutTraj(isnan(OutTraj)) = 0;
        clear OutTraj_x OutTraj_y
    end
end

% Define Exponent of exp
if(~Ift_flag)          
    Expy = -2*pi*1i;
else
    Expy = 2*pi*1i;
end

% Define Output Size
NOut = size(OutTraj,2);


%% 1. FT


sft1_Oper = zeros([size(OutTraj,2) size(InTraj,2)]);
for j=1:NOut
    sft1_Oper(j,:) = exp(Expy*( OutTraj(j)*InTraj ));
end



%% Normalize



if(Ift_flag)
   sft1_Oper = sft1_Oper / (size(InTraj,2)); 
end

