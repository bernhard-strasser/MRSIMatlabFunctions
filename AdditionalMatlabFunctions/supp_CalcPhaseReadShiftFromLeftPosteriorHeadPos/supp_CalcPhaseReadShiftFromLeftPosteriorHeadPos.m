function [Pos_PRS,Pars] = supp_CalcPhaseReadShiftFromLeftPosteriorHeadPos(Pars)
%
% ZerofillOrCutkSpace Zerofill or Cut k-Space Data
%
% This function was written by Bernhard Strasser, March 2019.
%
%
% The function either zerofills or cuts data in k-space. For zerofilling time data 
% (zerofilling only one-sided instead of on both sides like in k-space), use "Zerofilling_Spectral".
%
%
% OutArray = ZerofillOrCutkSpace(OutArray,Zerofill_To,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied.
% -     Zerofill_To                  ...    Array to which data should be zerofilled or cut. E.g. size(OutArray) = [32 64 64 512], Zerofill_To = [32 128 128 512]. 
% -     PerformFFT_flag              ...    If it is 1, the image gets Fourier transformed to k-space before applying the filter, 
%                                           and transformed back to image domain afterwards
%
% Output:
% -     OutArray                     ...    The filtered/masked output array
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%     Pars = read_ascconv('/autofs/space/carpathia_001/users/ovidiu/Data/7T/MeasuredData/MRSI_vol_7T_12Feb2019/ASE_TE78/MRSI_7T_VOL_12FEB2019.MR.INVESTIGATORS_OVIDIU.0013.0001.2019.02.12.10.59.47.554500.439749083.IMA');


% Compute direction cosine from x,y and z components of slice normal vector
[Pars.PhaseNormalVector, Pars.ReadNormalVector] = compute_dircos_1_3([Pars.SliceNormalVector_x(1) Pars.SliceNormalVector_y(1) Pars.SliceNormalVector_z(1)],Pars.InPlaneRotation);
Pars.SliceNormalVector = [Pars.SliceNormalVector_x(1) Pars.SliceNormalVector_y(1) Pars.SliceNormalVector_z(1)];    

MinusVec1 = [-1 -1 1];        % MinusVec are here only for "tuning" signs of the final rotation matrix (RotMat). The values of MinusVecs are 
MinusVec2 = [1 1 -1];         % empirical, thus there might exist a dataset, which will have signs of RotMat uncorrect. However, this setup
MinusVec3 = [-1 -1 +1];       % works for all tested datasets. The uncorrect signs of direction cosines might be caused by the fact, that the attribute
                              % Image Orientation (Patient) tag(0020,0037) is not used in computation of dircos. Instead the Siemens private Tag (0029, 1020) 
                              % is used

Pars.ReadNormalVector = Pars.ReadNormalVector .* MinusVec1;
Pars.PhaseNormalVector = Pars.PhaseNormalVector .* MinusVec2;
Pars.SliceNormalVector = Pars.SliceNormalVector .* MinusVec3;

% Create rotation matrix
RotMat = cat(1,Pars.ReadNormalVector,Pars.PhaseNormalVector,Pars.SliceNormalVector);

% Do we need to reverse x- and y-coordinates (taking the negative values)?
Pos_LPH = [Pars.Pos_Sag(1),Pars.Pos_Cor(1),Pars.Pos_Tra(1)];

Pos_PRS = RotMat * transpose(Pos_LPH); 

