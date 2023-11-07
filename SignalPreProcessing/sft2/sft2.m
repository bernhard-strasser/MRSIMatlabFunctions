function A_sft2 = sft2(A,InTraj,OutTraj,Ift_flag)
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
% A_sft2 = sft2(A,Ifft_flag)
%
% Input: 
% -         A                           ...     Input data.
% -         Ift_flag                    ...     Flag of 0 or 1 determining if an inverse Fourier transform or normal should be performed.
%
% Output:
% -         A_sft2                      ...     Fourier transformed result
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ?







%% 0. Preparations

% Input Sizes
N1 = size(A,1);
N2 = size(A,2);
N1_floor = floor(N1/2);
N2_floor = floor(N2/2);

if(~Ift_flag)          
    InScale = [N1 N2];                 % If doing an FT, assume InTraj is k-Space trajectory --> Scale k-space traj from [-0.5,0.5)
    OutScale = [1 1];                  % And OutTraj is i-Space trajectory --> Scale i-space traj from [-N/2, N/2)
    Expy = -2*pi*1i;
else
    InScale = [1 1];
    OutScale = [N1 N2];
    Expy = 2*pi*1i;
end

if(~exist('InTraj','var') || isempty(InTraj))
	[InTraj_x, InTraj_y] = ndgrid( ((1:N1) - N1_floor - 1)/InScale(1), ((1:N2) - N2_floor - 1)/InScale(2) );
	InTraj = cat(2,InTraj_x(:),InTraj_y(:));
    InTraj(isnan(InTraj)) = 0;
    clear InTraj_x InTraj_y
end
if(~exist('OutTraj','var') || isempty(OutTraj))
	[OutTraj_x, OutTraj_y] = ndgrid( ((1:N1) - N1_floor - 1)/OutScale(1), ((1:N2) - N2_floor - 1)/OutScale(2) );
	OutTraj = cat(2,OutTraj_x(:),OutTraj_y(:));
    OutTraj(isnan(OutTraj)) = 0;
    clear OutTraj_x OutTraj_y
end


% Define Output Size
NOut = size(OutTraj,1);

% Reshape Input
A = reshape(A,[size(InTraj,1) numel(A)/size(InTraj,1)]);



%% 1. FT

A_sft2 = zeros([size(OutTraj,1) size(A,2)]);
for j=1:NOut
    Expy2 = exp(Expy*( OutTraj(j,1)*InTraj(:,1) + OutTraj(j,2)*InTraj(:,2) ));
    Expy2 = repmat(Expy2,[1 size(A,2)]);
    A_sft2(j,:) = sum(   A.*Expy2  );     
end



%% Normalize


% A_sft2 = A_sft2 / (size(InTraj,1)); 

if(Ift_flag)
   A_sft2 = A_sft2 / (size(InTraj,1)); 
end

