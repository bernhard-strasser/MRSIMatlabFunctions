function A_sft2 = sft2_Cart(A,Ift_flag)
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



N1 = size(A,1);
N2 = size(A,2);



%% 1. FT

A_sft2 = zeros([N1 N2]); 

if(~Ift_flag)
    Expy = -2*pi*1i;
else
    Expy = 2*pi*1i;
end
    
N1_floor = floor(N1/2);
N1_ceil = ceil(N1/2);
N2_floor = floor(N2/2);
N2_ceil = ceil(N2/2);

% for j1=(-N1_floor):(N1_ceil-1)
%     for j2=(-N2_floor):(N2_ceil-1)
%         for k1=(-N1_floor):(N1_ceil-1)
%             for k2=(-N2_floor):(N2_ceil-1)
%                 A_sft2(j1+N1_floor+1,j2+N2_floor+1) = ...
%                 A_sft2(j1+N1_floor+1,j2+N2_floor+1) + A(k1+N1_floor+1,k2+N2_floor+1)*exp(Expy*((j1*k1)/N1+(j2*k2)/N2)); 
%             end
%         end
%     end
% end


for j1=(-N1_floor):(N1_ceil-1)
    for j2=(-N2_floor):(N2_ceil-1)
        k1=(-N1_floor):(N1_ceil-1);
        k2=(-N2_floor):(N2_ceil-1);
        Expy2 = exp(Expy*( myrepmat((j2*k2)/N2,[N1 N2],[0 2]) + myrepmat((j1*k1)/N1,[N1 N2],[2 0]) ));
        A_sft2(j1+N1_floor+1,j2+N2_floor+1) = sum(sum(   A.*Expy2  )); 
            
    end
end



if(Ift_flag)
   A_sft2 = A_sft2 / (N1*N2); 
end

