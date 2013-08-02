function [OutData,weights]=opengrappa_MRSI_1_2(InData, ACS, R, ApplyAlongDim) 
% 
% opengrappa_MRSI_x_y Reconstruct Undersampled MRSI Data
% 
% This code is based on the teaching version of the opengrappa function of Felix Breuer that was provided 
% on the parallel imaging workshop in Wuerzburg, 2012, which in turn is based on the opengrappa of Mark Griswold.
%
% 22.06.2008 Felix Breuer (breuer@mr-bavaria.de)
% November 2012 - January 2013 Bernhard Strasser
%
%
%
%  [OutData,weights] = opengrappa(InData,ACS,R);
%       
%   Input:      
% -     InData             Undersampled Data            (size: [#coils, k_x/R, k_y, slc, vecSize] or [#coils, k_x, k_y/R, slc, vecSize])
% -     ACS                AutoCalibration Signal       (size: [#coils, k_y, k_x, slices])     
% -     R                  Acceleration factor          (integer)
% -     ApplyAlongDim      Along which dimension        2 for along x-dimension, 3 for along y-dimension
%                          GRAPPA should be applied
%   Output:
% -     OutData            Reconstructed Output Data    (size: [#coils, k_x, k_y, slc, vecSize])   
% -     weights            Weights for Reconstruction              
%              
%
%   Some things to think about when using this code:
%
%           -The ACS lines used for reconstruction are NOT included in the final reconstructed
%           data sets. If you want to do this, you can do it after reconstruction. Please check
%           that they are lined up with the reconstructed data. I have not checked this.
%
%           -Since the ACS lines are not included, feel free to use a different imaging sequence
%           for acquisition of the ACS lines. We have seen some advantage in doing this when
%           the sequence used for the reduced acquisition is very flow sensitive, for example.
%           This can also be faster in many cases.
%
%           -The 4x5 block size is normally optimal for most applications, but if you need  
%           something faster, you should try the 4x3 code which uses a smaller block size.
%
%           -The matrix problems here are normally so overdetermined that numerical conditioning 
%           is normally not necessary. You could try something out on your problem and see what
%           you get. I have never seen any real advantage.
%
%           -The cyclic boundary condition will not work as programmed for radial and spiral
%           reconstructions. You will have to reprogram this yourself.
%           
%   Please read the license text at the bottom of this program. By using this program, you 
%   implicity agree with the license. 
%
%   The main points of the license:
%
%   1) This code is strictly for non-commercial applications. The code is protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
%



%% 0. Preparation

% Assign standard values to variables if nothing is passed to function.
if(~exist('ApplyAlongDim','var'))
    ApplyAlongDim = 2;
end 




%% 1. Apply Along ROW dimension

if(ApplyAlongDim == 2)

    %% 1.1 Preparations
    [nChannel,nx_red,ny, nSlice, nTime] = size(InData);  % Get the size of both the input data and the autocalibration data
    nx = nx_red * R;
    [nChannel_ACS,nx_ACS,ny_ACS, nSlice_ACS]=size(ACS);

    if nChannel_ACS~=nChannel
        disp('Error! The number of coils has to be the same for both inputs!')
        return;
    end

    fprintf('GRAPPA: \n')

    kernelsize_x = 4;                         % should be odd                  
    kernelsize_y = 5;                         % should be even 

    dy = floor(kernelsize_y/2);
    dx = (kernelsize_x/2-1)*R; 


    
    
    %% 1.2 Calculate weights

    tic; fprintf('Calculating weights')

    % Prepare source and target matrix
    % Number of source points -> nChannel*nsry*nsrcx
    % Number of target points -> nChannel * (R-1)
    % number of kernel repetitions in ACS data: (ny_ACS-(kernelsize_y-1)*R)*(nx_ACS-(kernelsize_x-1))

    SourcePoints_ACS = zeros(nChannel*kernelsize_x*kernelsize_y,(nx_ACS-(kernelsize_x-1)*R)*(ny_ACS-(kernelsize_y-1)));
    TargetPoints_ACS = zeros(nChannel*(R-1),(nx_ACS-(kernelsize_x-1)*R)*(ny_ACS-(kernelsize_y-1)));



    %   Simple example at R=2 and kernel size 4x5:               
    %
    %                   1 2 3 4 5
    %
    %                   0 0 0 0 0    1   
    %                   - - - - -    2
    %                   0 0 0 0 0    3       
    %                   - - X - -    4         
    %                   0 0 0 0 0    5
    %                   - - - - -    6
    %                   0 0 0 0 0    7
    %                   - - - - -    8  
    %
    %   The circles are the source points, and the X are the target points, - are undersampled k-points)
    %   Note: All source points in all coils are used to fit one target point in each individual coil 




    cnt = 0;                          %This is a very lazy counter. Could be done much faster. 

    for yind=dy+1:ny_ACS-dy,
        for xind=1:nx_ACS-(kernelsize_x-1)*R,
            cnt=cnt+1;

            % These are the source points (#coils*kernelsize_y*kernelsize_x) 
            SourcePoints_ACS(:,cnt) = reshape(ACS(:,xind:R:xind+(kernelsize_x-1)*R,yind-dy:yind+dy), nChannel*kernelsize_x*kernelsize_y,1);                   

            % these are the taget points (#coils*(R-1))
            TargetPoints_ACS(:,cnt) = reshape(ACS(:,xind+1+dx:xind+dx+R-1,yind), nChannel*(R-1),1);

        end
    end


    weights = TargetPoints_ACS * pinv(SourcePoints_ACS);  % find weights by fitting the source data to target data  


    fprintf('... %f sec \n',toc)                                                          

    

    %% 1.3 Apply weights

    tic; fprintf('Applyig weights ')


        
    Reco_dummy = zeros(nChannel,nx+2*dx+1,ny+2*dy,nSlice,nTime);                                                       % prepare matrix for convolution   
    Reco_dummy(:,dx+1:R:nx+dx,dy+1:end-dy,:,:) = InData;                                        % Write undersampled data into zero-padded matrix 

    for yind = dy+1:ny+dy, 
        for xind= 1:R:nx,
            SourcePoints_Reco=reshape(Reco_dummy(:,xind:R:xind + (kernelsize_x-1)*R,yind-dy:yind+dy,:,:), [nChannel*kernelsize_x*kernelsize_y nSlice*nTime]);
            Reco_dummy(:,xind+dx+1:xind+dx+R-1,yind,:,:)=reshape(weights*SourcePoints_Reco, [nChannel (R-1) nSlice nTime]);           %Apply weights to source points
        end
    end

    OutData = Reco_dummy(:,dx+1:nx+dx,dy+1:ny+dy,:,:);                                           %Crop out the good data.



    fprintf('... %f sec \n',toc)







%     tic; fprintf('Applyig weights ')
%     OutData = zeros([size(InData,1) size(InData,2)*R size(InData,3) size(InData,4) size(InData,5)]);
%         
%         Reco_dummy = zeros(nChannel,ny+2*dy+1,nx+2*dx,1, size(InData,5));                                                       % prepare matrix for convolution   
%         Reco_dummy(:,dy+1:R:ny+dy,dx+1:end-dx,1,:) = InData;                                        % Write undersampled data into zero-padded matrix 
% 
%         for xind = dx+1:nx+dx, 
%             for yind= 1:R:ny,
% 
%                 SourcePoints_Reco=reshape(Reco_dummy(:,yind:R:yind + (kernelsize_y-1)*R,xind-dx:xind+dx),nChannel*kernelsize_y*kernelsize_x,1);
% 
%                 Reco_dummy(:,yind+dy+1:yind+dy+R-1,xind)=reshape(weights*SourcePoints_Reco,[nChannel (R-1)]);           %Apply weights to source points
% 
%             end
%         end
%         OutData = Reco_dummy(:,dy+1:ny+dy,dx+1:nx+dx);                                           %Crop out the good data.
% 
% 
%     fprintf('... %f sec \n',toc)









%% 2. Apply Along COL Dimension

elseif(ApplyAlongDim == 3)
  
    
    %% 2.1 Preparations
    
    [nChannel,nx,ny_red, nSlice, nTime]=size(InData);
    ny = ny_red * R;
    
    [nChannel_ACS,nx_ACS,ny_ACS,nSlice_ACS]=size(ACS);     %Get the size of both the input data and the ref data

    if nChannel_ACS~=nChannel
        disp('Error! The number of coils has to be the same for both inputs!')
        return;
    end

    fprintf('GRAPPA: \n')

    kernelsize_x = 5;                         % should be odd                  
    kernelsize_y = 4;                         % should be even 

    dx = floor(kernelsize_x/2);
    dy = (kernelsize_y/2-1)*R; 

    

    %% 2.2   Calculate weights
   
    tic; fprintf('Calculating weights')

    % Prepare source and target matrix
    % Number of source points -> nChannel*nsrx*nsrcy
    % Number of target points -> nChannel * (R-1)
    % number of kernel repetitions in ACS data: (nx_ACS-(kernelsize_x-1)*R)*(ny_ACS-(kernelsize_y-1))
    %y = nChannel*kernelsize_x*kernelsize_y;
    %x = (nx_ACS-(kernelsize_x-1)*R)*(ny_ACS-(kernelsize_y-1));

    SourcePoints_ACS = zeros(nChannel*kernelsize_y*kernelsize_x,(ny_ACS-(kernelsize_y-1))*(nx_ACS-(kernelsize_x-1)*R));
    TargetPoints_ACS = zeros(nChannel*(R-1),(ny_ACS-(kernelsize_y-1))*(nx_ACS-(kernelsize_x-1)*R));

    
    
    
    %   Simple example at R=2 and kernel size 5x4:            
    %
    %                   1 2 3 4 5 6 7 8
    %
    %                   O - O - O - O -     1 
    %                   O - O - O - O -     2  
    %                   O - O X O - O -     3
    %                   O - O - O - O -     4
    %                   O - O - O - O -     5    
    %
    %   The circles are the source points, and the X are the target points, the - are undersampled data points 
    %   Note: All source points in all coils are used to fit one target point in each individual coil.

    
    
    
    cnt = 0;                          %This is a very lazy counter. Could be done much faster. 

    for yind=1:ny_ACS-(kernelsize_y-1)*R,
        for xind=dx+1:nx_ACS-dx,
            cnt=cnt+1;
            % These are the source points (#coils*kernelsize_y*kernelsize_x) 
            SourcePoints_ACS(:,cnt) = reshape(ACS(:,xind-dx:xind+dx, yind:R:yind+(kernelsize_y-1)*R), [nChannel*kernelsize_x*kernelsize_y 1]);    

            % these are the taget points (#coils*(R-1))
            TargetPoints_ACS(:,cnt) = reshape(ACS(:, xind, yind+1+dy:yind+dy+R-1), [nChannel*(R-1) 1] );

        end
    end

    weights=TargetPoints_ACS*pinv(SourcePoints_ACS);  % find weights by fitting the source data to target data
    clear SourcePoints_ACS TargetPoints_ACS
    
    fprintf('... %f sec \n',toc)                                                          

  
    
    
    %% 2.3   Apply weights


    tic; fprintf('Applyig weights ')

    Reco_dummy = zeros(nChannel,nx + 2*dx,ny + 2*dy + 1, nSlice, nTime);                                                             % prepare matrix for convolution     
    Reco_dummy(:,dx+1:end-dx,dy+1:R:ny+dy,:,:) = InData;                                                                             % Write undersampled data into zero-padded matrix 

    for yind = 1:R:ny, 
        for xind = dx+1:nx+dx
            SourcePoints_Reco=reshape(Reco_dummy(:,xind-dx:xind+dx, yind:R:yind + (kernelsize_y-1)*R,1,:), [nChannel*kernelsize_x*kernelsize_y nSlice*nTime]);
            Reco_dummy(:,xind,yind+dy+1:yind+dy+R-1,:,:)=reshape(weights*SourcePoints_Reco,[nChannel (R-1) nSlice nTime]);           %Apply weights to source points
        end
    end
    
    OutData = Reco_dummy(:,dx+1:nx+dx,dy+1:ny+dy,:,:);                                                                               % Crop out the good data.

    fprintf('... %f sec \n',toc)
    
    
    
    
else
    disp('Error! The the type parameter value is wrrong!')
    OutData = 0;
    weights = 0;
    return;
end








