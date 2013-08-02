function [recon,ws]=opengrappa_time(sig, acs, af, type) 

% This code is written by Felix Breuer based on opengrappa by Mark Griswold 
%
% 22.06.2008 Felix Breuer (breuer@mr-bavaria.de)
%
%
%  [recon,sigrecon,ws,ws_img,g] = opengrappa(sig,acs,af);
%
%   This is a teaching version of the GRAPPA reconstruction code.
%       
%   IN:         sig                reduced data set         (#coils, Ky./af, Kx)
%               acs                autocalibration lines    (#coils, #acs lines, Kx-acslines)    
%               af                 Acceleration factor      (integer)
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


if strcmp(type, 'row') == 1

    [nc,ny_red,nx, ns, nTime]=size(sig);
    [nc_acs,nyacs,nxacs]=size(acs);     %Get the size of both the input data and the ref data
    ny = ny_red *af;

    if nc_acs~=nc
        disp('Error! The number of coils has to be the same for both inputs!')
        return;
    end

    fprintf('GRAPPA: \n')

    srcx = 5;                         % should be odd                  
    srcy = 4;                         % should be even 


    dx = floor(srcx/2);
    dy = (srcy/2-1)*af; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Calculate weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic; fprintf('Calculating weights')

    % Prepare source and target matrix
    % Number of source points -> nc*nsry*nsrcx
    % Number of target points -> nc * (af-1)
    % number of kernel repetitions in ACS data: (nyacs-(srcy-1)*af)*(nxacs-(srcx-1))

    src = zeros(nc*srcy*srcx,(nyacs-(srcy-1)*af)*(nxacs-(srcx-1)));
    trg = zeros(nc*(af-1),(nyacs-(srcy-1)*af)*(nxacs-(srcx-1)));



    %   Simple example at af=2 and kernel size 5x4:            
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
        for yind=1:nyacs-(srcy-1)*af,
            cnt=cnt+1;

            % These are the source points (#coils*srcy*srcx) 
            src(:,cnt) = reshape(acs(:,yind:af:yind+(srcy-1)*af,xind-dx:xind+dx),nc*srcy*srcx,1);                   

            % these are the taget points (#coils*(af-1))
            trg(:,cnt) = reshape(acs(:,yind+1+dy:yind+dy+af-1,xind),nc*(af-1),1);

        end
    end




    %ws=trg*pinv(src);  % find weights by fitting the source data to target data  
    ws=trg*pinv(squeeze(src(:,:,:, 1, 1)));

    t=toc; fprintf('... %f sec \n',t)                                                          

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Apply weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic; fprintf('Applyig weights ')
    for i=1:nTime,
        sigrecon = zeros(nc,ny+2*dy+1,nx+2*dx);                                                          % prepare matrix for convolution   

        sigrecon(:,dy+1:af:ny+dy,dx+1:end-dx) = sig(:,:,:, 1, i);                                        % Write undersampled data into zero-padded matrix 

        for xind = dx+1:nx+dx, 
            for yind= 1:af:ny,

                src=reshape(sigrecon(:,yind:af:yind + (srcy-1)*af,xind-dx:xind+dx),nc*srcy*srcx,1);

                sigrecon(:,yind+dy+1:yind+dy+af-1,xind)=reshape(ws*src,[nc (af-1)]);           %Apply weights to source points

            end
        end
        sigrecon = sigrecon(:,dy+1:ny+dy,dx+1:nx+dx);                                           %Crop out the good data.

        recon_tmp=fftshift(fftshift(fft(fft(ifftshift(ifftshift(sigrecon,2),3),[],2),[],3),2),3);
        recon(:, :, :, 1, i) = flipdim(recon_tmp,2);
    end

    t=toc;
    fprintf('... %f sec \n',t)
else
if strcmp(type, 'col') == 1
  
    [nc,ny,nx_red, ns, nTime]=size(sig);
    [nc_acs,nyacs,nxacs]=size(acs);     %Get the size of both the input data and the ref data
    nx = nx_red *af;

    if nc_acs~=nc
        disp('Error! The number of coils has to be the same for both inputs!')
        return;
    end

    fprintf('GRAPPA: \n')

    srcx = 4;                         % should be odd                  
    srcy = 5;                         % should be even 

    dy = floor(srcy/2);
    dx = (srcx/2-1)*af; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Calculate weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic; fprintf('Calculating weights')

    % Prepare source and target matrix
    % Number of source points -> nc*nsry*nsrcx
    % Number of target points -> nc * (af-1)
    % number of kernel repetitions in ACS data: (nyacs-(srcy-1)*af)*(nxacs-(srcx-1))
    %x = nc*srcy*srcx;
    %y = (nyacs-(srcy-1)*af)*(nxacs-(srcx-1));

    src = zeros(nc*srcy*srcx,(nyacs-(srcy-1))*(nxacs-(srcx-1)*af));
    trg = zeros(nc*(af-1),(nyacs-(srcy-1))*(nxacs-(srcx-1)*af));

    %   Simple example at af=2 and kernel size 5x4:            
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

    for xind=1:nxacs-(srcx-1)*af,
        for yind=dy+1:nyacs-dy,
            cnt=cnt+1;
            % These are the source points (#coils*srcy*srcx) 
            src(:,cnt) = reshape(acs(:,yind-dy:yind+dy, xind:af:xind+(srcx-1)*af), nc*srcy*srcx, 1);    
            %src(:,cnt) = reshape(acs(:,yind:af:yind+(srcy-1)*af,xind-dx:xind+dx),nc*srcy*srcx,1);               

            % these are the taget points (#coils*(af-1))
            trg(:,cnt) = reshape(acs(:, yind, xind+1+dx:xind+dx+af-1),nc*(af-1),1);
            %trg(:,cnt) = reshape(acs(:,yind+1+dy:yind+dy+af-1,xind),nc*(af-1),1);

        end
    end

    %ws=trg*pinv(src);  % find weights by fitting the source data to target data  
    ws=trg*pinv(squeeze(src(:,:,:, 1, 1)));

    t=toc; fprintf('... %f sec \n',t)                                                          

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Apply weights
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic; fprintf('Applyig weights ')
    for i = 1:nTime,
        sigrecon = zeros(nc,ny+2*dy,nx+2*dx+ 1);                                                         % prepare matrix for convolution     
        sigrecon(:,dy+1:end-dy,dx+1:af:nx+dx) = sig(:,:,:, 1, i);                                        % Write undersampled data into zero-padded matrix 
        %sigrecon(:,dy+1:af:ny+dy,dx+1:end-dx) = sig(:,:,:, 1, i);

        for xind = 1:af:nx, 
            for yind = dy+1:ny+dy
                src=reshape(sigrecon(:,yind-dy:yind+dy, xind:af:xind + (srcx-1)*af),nc*srcy*srcx,1);
                %src=reshape(sigrecon(:,yind:af:yind + (srcy-1)*af,xind-dx:xind+dx),nc*srcy*srcx,1);
                sigrecon(:,yind,xind+dx+1:xind+dx+af-1)=reshape(ws*src,[nc (af-1)]);           %Apply weights to source points
                %sigrecon(:,yind+dy+1:yind+dy+af-1,xind)=reshape(ws*src,[nc (af-1)]);
            end
        end
        sigrecon = sigrecon(:,dy+1:ny+dy,dx+1:nx+dx);                                           %Crop out the good data.
        recon_tmp=fftshift(fftshift(fft(fft(ifftshift(ifftshift(sigrecon,2),3),[],2),[],3),2),3);
        recon(:, :, :, 1, i) = flipdim(recon_tmp,2);
    end

    t=toc;
    fprintf('... %f sec \n',t)
else
    disp('Error! The the type parameter value is wrrong!')
    recon = 0;
    ws = 0;
    return;
end
end