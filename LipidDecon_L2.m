function OutData = LipidDecon_L2(InData,param)
%-----------------------------------------------------------------------
% Based on M. Lustig's Conjugate Gradient minimizer
%
% given k-space measurements d, and a fourier operator F the function 
% finds the image m that minimizes:
%
% Phi(x) = || F * m - k_space ||^2 +       
%          lambda * sum_over((x,y) in brain_mask){ ||L * m(x,y)||_1 }
%
% the optimization method used is non linear conjugate gradient with 
% backtracking line-search.
% 
% params.FT : DFT operator 
%
% params.data :  k-space data of lipid contaminated image
%
% params.Bmask : brain mask
% params.Lipid : lipid basis functions
%
%-------------------------------------------------------------------------
	Lipid_inv = inv( eye(size(InData,3)) + param.beta * (param.Lipid * param.Lipid') );
	OutData = InData;

	for ay = 1:size(OutData,1)
		for cey = 1:size(OutData,2)

			if (param.Bmask(ay,cey))

				xi = InData(ay,cey,:);
				xi_L2 = Lipid_inv * xi(:);
				OutData(ay,cey,:) = xi_L2;

			end

		end
	end

	OutData = OutData .* myrepmat(param.Bmask,size(OutData));

return;