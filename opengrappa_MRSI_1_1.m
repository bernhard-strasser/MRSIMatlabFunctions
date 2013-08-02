function [recon,ws]=opengrappa_MRSI_1_1(sig, acs, R, type) 

% This code is written by Felix Breuer based on opengrappa by Mark Griswold 
%
% 22.06.2008 Felix Breuer (breuer@mr-bavaria.de)
%
%
%  [recon,sigrecon,ws,ws_img,g] = opengrappa(sig,acs,R);
%
%   This is a teaching version of the GRAPPA reconstruction code.
%       
%   IN:         sig                reduced data set         (#coils, Ky./R, Kx)
%               acs                autocalibration lines    (#coils, #acs lines, Kx-acslines)    
%               R                  Acceleration factor      (integer)
%
%   OUT:        recon              Reconstructed images     (# coils, Ny, Nx)   
%               ws                 Weight sets              
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
%           -Finally, the coil index is the first index simply for historical reasons. I will
%           eventually reprogram this. Use the 'permute' command to reorder your matrices if you
%           have the coil index last. 
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














if(strcmp(type, 'row') == 1)

    [nc,ny_red,nx, ns, nTime]=size(sig);
    ny = ny_red * R;
    [nc_acs,nyacs,nxacs]=size(acs);     %Get the size of both the input data and the ref data

    if nc_acs~=nc
        disp('Error! The number of coils has to be the same for both inputs!')
        return;
    end

    fprintf('GRAPPA: \n')

    srcx = 5;                         % should be odd                  
    srcy = 4;                         % should be even 


    dx = floor(srcx/2);
    dy = (srcy/2-1)*R; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Calculate weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic; fprintf('Calculating weights')

    % Prepare source and target matrix
    % Number of source points -> nc*nsry*nsrcx
    % Number of target points -> nc * (R-1)
    % number of kernel repetitions in ACS data: (nyacs-(srcy-1)*R)*(nxacs-(srcx-1))

    src = zeros(nc*srcy*srcx,(nyacs-(srcy-1)*R)*(nxacs-(srcx-1)));
    trg = zeros(nc*(R-1),(nyacs-(srcy-1)*R)*(nxacs-(srcx-1)));



    %   Simple example at R=2 and kernel size 5x4:            
    %
    %                   -O-O-O-O-O-     1   
    %                   - - - - - -     2
    %                   -O-O-O-O-O-     3       
    %                   - - -X- - -     4         
    %                   -O-O-O-O-O-     5
    %                   - - - - - -     6
    %                   -O-O-O-O-O-     7
    %                   - - - - - -     8
    %
    %   The circles are the source points, and the X are the target points.
    %   Note: All source points in all coils are used to fit one target point in each individual coil 




    cnt = 0;                          %This is a very lazy counter. Could be done much faster. 

    for xind=dx+1:nxacs-dx,
        for yind=1:nyacs-(srcy-1)*R,
            cnt=cnt+1;

            % These are the source points (#coils*srcy*srcx) 
            src(:,cnt) = reshape(acs(:,yind:R:yind+(srcy-1)*R,xind-dx:xind+dx),nc*srcy*srcx,1);                   

            % these are the taget points (#coils*(R-1))
            trg(:,cnt) = reshape(acs(:,yind+1+dy:yind+dy+R-1,xind),nc*(R-1),1);

        end
    end




    %ws=trg*pinv(src);  % find weights by fitting the source data to target data  
    ws=trg*pinv(squeeze(src(:,:,:, 1, 1)));

    fprintf('... %f sec \n',toc)                                                          

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Apply weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic; fprintf('Applyig weights ')
    recon = zeros([size(sig,1) size(sig,2)*R size(sig,3) size(sig,4) size(sig,5)]);
    for t=1:nTime
        
        sigrecon = zeros(nc,ny+2*dy+1,nx+2*dx);                                                       % prepare matrix for convolution   
        sigrecon(:,dy+1:R:ny+dy,dx+1:end-dx) = sig(:,:,:, 1, t);                                        % Write undersampled data into zero-padded matrix 

        for xind = dx+1:nx+dx, 
            for yind= 1:R:ny,

                src=reshape(sigrecon(:,yind:R:yind + (srcy-1)*R,xind-dx:xind+dx),nc*srcy*srcx,1);

                sigrecon(:,yind+dy+1:yind+dy+R-1,xind)=reshape(ws*src,[nc (R-1)]);           %Apply weights to source points

            end
        end
        recon(:,:,:,:,t) = sigrecon(:,dy+1:ny+dy,dx+1:nx+dx);                                           %Crop out the good data.

    end

    fprintf('... %f sec \n',toc)







%     tic; fprintf('Applyig weights ')
%     recon = zeros([size(sig,1) size(sig,2)*R size(sig,3) size(sig,4) size(sig,5)]);
%         
%         sigrecon = zeros(nc,ny+2*dy+1,nx+2*dx,1, size(sig,5));                                                       % prepare matrix for convolution   
%         sigrecon(:,dy+1:R:ny+dy,dx+1:end-dx,1,:) = sig;                                        % Write undersampled data into zero-padded matrix 
% 
%         for xind = dx+1:nx+dx, 
%             for yind= 1:R:ny,
% 
%                 src=reshape(sigrecon(:,yind:R:yind + (srcy-1)*R,xind-dx:xind+dx),nc*srcy*srcx,1);
% 
%                 sigrecon(:,yind+dy+1:yind+dy+R-1,xind)=reshape(ws*src,[nc (R-1)]);           %Apply weights to source points
% 
%             end
%         end
%         recon = sigrecon(:,dy+1:ny+dy,dx+1:nx+dx);                                           %Crop out the good data.
% 
% 
%     fprintf('... %f sec \n',toc)












elseif(strcmp(type, 'col') == 1)
  
    [nc,ny,nx_red, ns, nTime]=size(sig);
    nx = nx_red * R;
    
    [nc_acs,nyacs,nxacs]=size(acs);     %Get the size of both the input data and the ref data

    if nc_acs~=nc
        disp('Error! The number of coils has to be the same for both inputs!')
        return;
    end

    fprintf('GRAPPA: \n')

    srcx = 4;                         % should be odd                  
    srcy = 5;                         % should be even 

    dy = floor(srcy/2);
    dx = (srcx/2-1)*R; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Calculate weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic; fprintf('Calculating weights')

    % Prepare source and target matrix
    % Number of source points -> nc*nsry*nsrcx
    % Number of target points -> nc * (R-1)
    % number of kernel repetitions in ACS data: (nyacs-(srcy-1)*R)*(nxacs-(srcx-1))
    %x = nc*srcy*srcx;
    %y = (nyacs-(srcy-1)*R)*(nxacs-(srcx-1));

    src = zeros(nc*srcy*srcx,(nyacs-(srcy-1))*(nxacs-(srcx-1)*R));
    trg = zeros(nc*(R-1),(nyacs-(srcy-1))*(nxacs-(srcx-1)*R));

    %   Simple example at R=2 and kernel size 5x4:            
    %
    %                   -O-O-O-O-O-     1   
    %                   - - - - - -     2
    %                   -O-O-O-O-O-     3       
    %                   - - -X- - -     4         
    %                   -O-O-O-O-O-     5
    %                   - - - - - -     6
    %                   -O-O-O-O-O-     7
    %                   - - - - - -     8
    %
    %   The circles are the source points, and the X are the target points.
    %   Note: All source points in all coils are used to fit one target point in each individual coil 

    cnt = 0;                          %This is a very lazy counter. Could be done much faster. 

    for xind=1:nxacs-(srcx-1)*R,
        for yind=dy+1:nyacs-dy,
            cnt=cnt+1;
            % These are the source points (#coils*srcy*srcx) 
            src(:,cnt) = reshape(acs(:,yind-dy:yind+dy, xind:R:xind+(srcx-1)*R), nc*srcy*srcx, 1);    
            %src(:,cnt) = reshape(acs(:,yind:R:yind+(srcy-1)*R,xind-dx:xind+dx),nc*srcy*srcx,1);               

            % these are the taget points (#coils*(R-1))
            trg(:,cnt) = reshape(acs(:, yind, xind+1+dx:xind+dx+R-1),nc*(R-1),1);
            %trg(:,cnt) = reshape(acs(:,yind+1+dy:yind+dy+R-1,xind),nc*(R-1),1);

        end
    end

    %ws=trg*pinv(src);  % find weights by fitting the source data to target data  
    ws=trg*pinv(squeeze(src(:,:,:, 1, 1)));

    fprintf('... %f sec \n',toc)                                                          

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Apply weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic; fprintf('Applyig weights ')
    recon = zeros([size(sig,1) size(sig,2) size(sig,3)*R size(sig,4) size(sig,5)]);
    for t = 1:nTime,
        sigrecon = zeros(nc,ny+2*dy,nx+2*dx+ 1);                                                         % prepare matrix for convolution     
        sigrecon(:,dy+1:end-dy,dx+1:R:nx+dx) = sig(:,:,:, 1, t);                                           % Write undersampled data into zero-padded matrix 
        %sigrecon(:,dy+1:R:ny+dy,dx+1:end-dx) = sig(:,:,:, 1, t);

        for xind = 1:R:nx, 
            for yind = dy+1:ny+dy
                src=reshape(sigrecon(:,yind-dy:yind+dy, xind:R:xind + (srcx-1)*R),nc*srcy*srcx,1);
                %src=reshape(sigrecon(:,yind:R:yind + (srcy-1)*R,xind-dx:xind+dx),nc*srcy*srcx,1);
                sigrecon(:,yind,xind+dx+1:xind+dx+R-1)=reshape(ws*src,[nc (R-1)]);           %Apply weights to source points
                %sigrecon(:,yind+dy+1:yind+dy+R-1,xind)=reshape(ws*src,[nc (R-1)]);
            end
        end
        recon(:,:,:,:,t) = sigrecon(:,dy+1:ny+dy,dx+1:nx+dx);                                           %Crop out the good data.

    end

    fprintf('... %f sec \n',toc)
    
    
    
    
else
    disp('Error! The the type parameter value is wrrong!')
    recon = 0;
    ws = 0;
    return;
end


