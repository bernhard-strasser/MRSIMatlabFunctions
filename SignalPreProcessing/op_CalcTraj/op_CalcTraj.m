function k = op_CalcTraj(MRStruct, Settings, Info)
    %
    % op_ReadAndRecoBorjanSpiralMRStruct.Data Read and reconstruct MRStruct.Data from Borjan Gagoski's Spiral MRSI Sequence
    %
    % This function was written by Bernhard Strasser, June 2019.
    %
    %
    % The function can read in Spiral MRSI MRStruct.Datas in the Siemens raw file format ".DAT" and performs
    % the reconstruction of the MRStruct.Data (Non-Uniform Slow FourierTransform etc.)
    %
    %
    % [MRStruct, Info] = read_csi_dat(file, DesiredSize,ReadInMRStruct.DataSets)
    %
    % Input: 
    % -         MRStruct          ...  An MRStruct   
    % -         TrajectoryFile    ...  The trajectory file containing information about the k-space spiral trajectory.
    % -         Settings                ...  Struct with fields
    %                                           Debug: The Debug-settings. Subfields:
    %                                                           ShowTrajs:  If you want to plot the spiral and Cartesian trajectories to get a feeling how you measured...
    %                                           io_ReadSpiralPars:  Settings for the function "io_ReadSpiralPars", see io_ReadSpiralPars.m
    %                                           ReadInTraj:         Settings for reading the trajectory, see io_ReadSpiralTraj.m.
    %                                           CalcOutTraj:        Settings for calculating the Cartesian trajectory, see sim_CalcCartTraj.m.
    %                                           NonCartReco:        Settings for the non-Cartesian MRSI Reco, see op_ReconstructNonCartMRMRStruct.Data.m.
    %
    % Output:
    % -         MRStruct            ...  Output-struct with fields
    %                                           RecoSteps:  All the settings of the performed reconstruction steps
    %                                           Par:        The original parameters of the read in MRStruct.Data
    %                                           RecoPar:    The parameters after reconstruction
    %                                           MRStruct.Data:       The reconstructed MRStruct.Data
    %                                           NoiseMRStruct.Data:  If available, noise with the same scale and size as the MRStruct.Data is produced. Useful for SNR calculations.
    %                                           OutTraj:    The Cartesian trajectory which would correspond to the Fourier transform of the Cartesian output image
    %                                                       (if we did k-space gridding, those would be the k-space points we would grid to).
    %                                           InTraj:     The (spiral) trajectory with which the MRStruct.Data were measured.
    %                                           
    % -         AdditionalOut           ...  Struct with additional output, e.g. Fourier Transform Operator, Density Compensation Function, etc.
    %
    %
    % Feel free to change/reuse/copy the function. 
    % If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
    % Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
    % File dependancy: ???

    % Further remarks:



    %% 0. Preparations




    if(~exist('Settings','var'))
       Settings = struct; 
    end

    
    %%
    
    
    tol= pi*1.0;

    gx_pos_left=unwrap(angle(MRStruct.Data(:,:,:,1,:,:,:,:)),tol,1);
    gx_pos_right=unwrap(angle(MRStruct.Data(:,:,:,5,:,:,:,:)),tol,1);
    gx_neg_left=unwrap(angle(MRStruct.Data(:,:,:,2,:,:,:,:)),tol,1);
    gx_neg_right=unwrap(angle(MRStruct.Data(:,:,:,6,:,:,:,:)),tol,1);
    gy_pos_left=unwrap(angle(MRStruct.Data(:,:,:,3,:,:,:,:)),tol,1);
    gy_pos_right=unwrap(angle(MRStruct.Data(:,:,:,7,:,:,:,:)),tol,1);
    gy_neg_left=unwrap(angle(MRStruct.Data(:,:,:,4,:,:,:,:)),tol,1);
    gy_neg_right=unwrap(angle(MRStruct.Data(:,:,:,8,:,:,:,:)),tol,1);
    % gz_pos_left=unwrap(angle(MRStruct.Data(:,:,:,5)),tol,1);
    % gz_pos_right=unwrap(angle(MRStruct.Data(:,:,:,11)),tol,1);
    % gz_neg_left=unwrap(angle(MRStruct.Data(:,:,:,6)),tol,1);
    % gz_neg_right=unwrap(angle(MRStruct.Data(:,:,:,12)),tol,1);
    % gx_pos_left=rescale(gx_pos_left,-pi,pi);

    % MRStruct.Data=permute(MRStruct.Data(:,:,1:2,1,:),[1 2 3 4 5]);
    % MRStruct.Data=reshape(MRStruct.Data,[size(MRStruct.Data,1)*size(MRStruct.Data,2)*size(MRStruct.Data,4)*size(MRStruct.Data,3) size(MRStruct.Data,5)]);

    offcentershift=MRStruct.Par.Pos_Tra; %from isocenter, see sourcecode
    Ax=gx_pos_left-gx_neg_left; Bx=gx_neg_right - gx_pos_right;
    Ay=gy_pos_left-gy_neg_left; By=gy_neg_right - gy_pos_right;
    % Az=gz_pos_left-gz_neg_left; Bz=gz_neg_right - gz_pos_right;

    k.x= squeeze_single_dim(1/4/offcentershift*(Ax+Bx),4,1);
    k.y= squeeze_single_dim(1/4/offcentershift*(Ay+By),4,1);





    
end



