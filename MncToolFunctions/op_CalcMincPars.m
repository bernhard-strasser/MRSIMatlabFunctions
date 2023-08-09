function MRStruct = op_CalcMincPars(MRStruct)
% read_MncFiles Read a simple binary 'raw' file, e.g. created by mnc2raw
%
% This function was written by Bernhard Strasser, September 2020.
%
%
% The function 
%
%
% image = read_MncFiles(MincFile,ROW,COL,SLC,precision)
%
% Input: 
% -        MRStruct                        ...     MR Structure with Par or RecoPar field.
%
%
% Output:
% -        MRStruct                        ...     Output will have added field MncPars.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None








%% 0. Preparations


if(isfield(MRStruct,'RecoPar'))
    InPars = MRStruct.RecoPar;
elseif(isfield(MRStruct,'Par'))
    InPars = MRStruct.Par;    
else
    error('MRStruct has neither subfield ''Pars'' nor ''RecoPars''. Cannot calculate MincPars.')
end


%%


	% Stepsize (Voxel Size)
	MRStruct.MncPars.StepRead = -InPars.FoV_Read(1) / InPars.nFreqEnc;		% Coordinate system is reversed in minc with respect to DICOM
	MRStruct.MncPars.StepPhase = -InPars.FoV_Phase(1) / InPars.nPhasEnc;    
	MRStruct.MncPars.StepSlice = InPars.FoV_Partition(1) / (InPars.nPartEnc * InPars.nSLC); 
	
	
	% Rename Position Fields
	[POS_X] = InPars.Pos_Sag;
	[POS_Y] = InPars.Pos_Cor;
	[POS_Z] = InPars.Pos_Tra;

    
	% Compute direction cosine from x,y and z components of slice normal vector
	[MRStruct.MncPars.PhaseNormalVector, MRStruct.MncPars.ReadNormalVector] = compute_dircos([InPars.SliceNormalVector_x(1) InPars.SliceNormalVector_y(1) InPars.SliceNormalVector_z(1)],InPars.InPlaneRotation);
	MRStruct.MncPars.SliceNormalVector = [InPars.SliceNormalVector_x(1) InPars.SliceNormalVector_y(1) InPars.SliceNormalVector_z(1)];    

	MinusVec1 = [-1 -1 1];        % MinusVec are here only for "tuning" signs of the final rotation matrix (RotMat). The values of MinusVecs are 
	MinusVec2 = [1 1 -1];         % empirical, thus there might exist a dataset, which will have signs of RotMat uncorrect. However, this setup
	MinusVec3 = [-1 -1 +1];       % works for all tested datasets. The uncorrect signs of direction cosines might be caused by the fact, that the attribute
								  % Image Orientation (Patient) tag(0020,0037) is not used in computation of dircos. Instead the Siemens private Tag (0029, 1020) 
								  % is used


	MRStruct.MncPars.ReadNormalVector = MRStruct.MncPars.ReadNormalVector .* MinusVec1;
	MRStruct.MncPars.PhaseNormalVector = MRStruct.MncPars.PhaseNormalVector .* MinusVec2;
	MRStruct.MncPars.SliceNormalVector = MRStruct.MncPars.SliceNormalVector .* MinusVec3;

	% Create rotation matrix
	RotMat = cat(1,MRStruct.MncPars.ReadNormalVector,MRStruct.MncPars.PhaseNormalVector,MRStruct.MncPars.SliceNormalVector);

	% Reverse x- and y- coordinates due to the reversed coordinate system of minc with respect to DICOM
	Pos = [-POS_X(1),-POS_Y(1),POS_Z(1)];

	% Convert position from DICOM world coordinates to MINC start values 
	% This approach is probably prone to extreme rotation of FOV and is not
	% universal, however for the tested datasets it yielded correct results
	Pos_Minc = RotMat * transpose(Pos); 
	% The following line is old, and I think wrong. What I think happened is: Michal fixed the shift
	% in the z-direction by half a voxel with the line below (effectively this line calculates
	% Pos_z - FoV_z/2 + FoV_z/(2*N_z). Then I figured out that the x- and y-positions have to be
	% shifted by half a voxel, and thought by analogy also the z-dimension has to be shifted, not
	% knowing that Michal did that already with the FoVHalf. Now it should be fixed:
	% The FoVHalf is defined "normally" also for z, and the half-voxel shift is done in the 3D-case.
	% For checks: See git commits #1383, #004f, #a9be
% 	FoVHalf = [InPars.FoV_Read(1)/2 InPars.FoV_Phase(1)/2 -InPars.FoV_Partition(1)/InPars.nPartEnc*(InPars.nPartEnc-1)/2]; 
    FoVHalf = [InPars.FoV_Read(1)/2 InPars.FoV_Phase(1)/2 -InPars.FoV_Partition(1)/2];  
	Pos_Minc = transpose(Pos_Minc) + FoVHalf;

	% Get from Center of Voxel (DICOM) to corner of voxel (minc) by subtracting half the voxel 
	MRStruct.MncPars.POS_X_FirstVoxel = Pos_Minc(1) + MRStruct.MncPars.StepRead/2;         % Be aware that StepRead and StepPhase are reversed and thus the sum is effectively a subtraction.
	MRStruct.MncPars.POS_Y_FirstVoxel = Pos_Minc(2) + MRStruct.MncPars.StepPhase/2;
	MRStruct.MncPars.POS_Z_FirstVoxel = Pos_Minc(3);

	
	if(InPars.ThreeD_flag)
		MRStruct.MncPars.POS_Z_FirstVoxel = MRStruct.MncPars.POS_Z_FirstVoxel + MRStruct.MncPars.StepSlice/2;     
	end




