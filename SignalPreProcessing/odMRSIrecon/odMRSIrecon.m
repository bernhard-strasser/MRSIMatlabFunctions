function csi_sum=odMRSIrecon(csi,B0map,t)
% computing only for one slice data, need to modify for 3D odMRSI


%% Housekeeping

% time vector
% nS  = Par.CSI.vecSize;
% dwt = Par.CSI.Dwelltimes(1);
% t   = (0:nS-1)*dwt/10^9;

S_map = size(B0map);
if(numel(S_map) < 3)
    S_map = [S_map 1];
end
S_csi = size(csi);


%% Interpolate B0-map so that its size is multiple integer of csi data
if(sum(mod(S_map,S_csi(1:3))) > 0)
	fac = round(S_map./S_csi(1:3)); 
    ny=fac(1)*S_csi(1);nx=fac(2)*S_csi(2);nz=fac(3)*S_csi(3); %% desired output dimensions
    [y, x, z] = ndgrid(linspace(1,size(B0map,1),ny), linspace(1,size(B0map,2),nx), linspace(1,size(B0map,3),nz));
    if(size(y,3) > 1)
        B0map=interp3(B0map,x,y,z,'cubic');
    else
        B0map=interp2(B0map,x,y,'cubic');
    end
    clear ny nx nz x y z
    S_map = size(B0map);
end


fac = S_map./S_csi(1:3);

% %% Denoising process (does not work well since it was done in nonlocal method)
% 
% for j = 1:S_csi(4)
%     for k = 1:S_csi(3)
%         csi_real0 = imgaussfilt(real(csi(:,:,k,j)));
%         csi_imag0 = imgaussfilt(imag(csi(:,:,k,j)));
%         csi(:,:,k,j) = csi_real0 + 1i*csi_imag0;
%     end
% end

%% Zerofill CSI Data
if(sum(fac) > numel(fac))
    csi = ZerofillOrCutkSpace(csi,S_csi.*[fac 1],1);
end



%% Frequency Align Data
csi = csi .* exp(myrepmat(2*pi*1i*B0map,size(csi)) .* myrepmat(t,size(csi)));



% h1 = fspecial3('gaussian', [3,3,3]);

% for i = 1:512
%     csi(:,:,:,i) = convn(csi(:,:,:,i), h1, 'same');
% end











% %% Boxcar with Gaussian Averaging (Need to adapt to more general case!)
% h1 = fspecial3('gaussian', [3,3,5]);
% csi_new = zeros([138,138,41,512]);
% csi_new(:,:,1:40,:) = csi;
% csi_sum = zeros(S_csi);
% for i = 1:512
%     for ro = 1:S_csi(1)
%         for co = 1:S_csi(2)
%             for zz = 1:S_csi(3)
%                 DP = csi_new((3*ro-2):3*ro,(3*co-2):3*co,(4*zz-3):(4*zz+1),i).*h1;
%                 csi_sum(ro,co,zz,i) = sum(DP(:));
%             end
%         end
%     end
%     
% end

%% General cases
fac_new = zeros(1,3);
S_map_new = zeros(size(S_map));
for ii = 1:3
    if mod(fac(ii),2) == 0
        fac_new(ii) = fac(ii) + 1;
        S_map_new(ii) = S_map(ii) + 1;
    else
        fac_new(ii) = fac(ii);
        S_map_new(ii) = S_map(ii);
    end
end

h1 = fspecial3('gaussian', fac_new);
csi_new = zeros([S_map_new,size(csi,4)]);
csi_new(1:size(csi,1),1:size(csi,2),1:size(csi,3),1:size(csi,4)) = csi;
csi_sum = zeros(S_csi);
if mod(fac(1),2) ~= 0 && mod(fac(2),2) ~= 0 && mod(fac(3),2) ~= 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):fac(1)*ro,(fac(2)*co-(fac(2)-1)):fac(2)*co,(fac(3)*zz-(fac(3)-1)):fac(3)*zz, i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
    
elseif mod(fac(1),2) ~= 0 && mod(fac(2),2) ~= 0 && mod(fac(3),2) == 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):fac(1)*ro,(fac(2)*co-(fac(2)-1)):fac(2)*co,(fac(3)*zz-(fac(3)-1)):(fac(3)*zz+1), i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
    
elseif mod(fac(1),2) ~= 0 && mod(fac(2),2) == 0 && mod(fac(3),2) ~= 0    
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):fac(1)*ro,(fac(2)*co-(fac(2)-1)):(fac(2)*co+1),(fac(3)*zz-(fac(3)-1)):fac(3)*zz, i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
    
elseif mod(fac(1),2) == 0 && mod(fac(2),2) ~= 0 && mod(fac(3),2) ~= 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):(fac(1)*ro+1),(fac(2)*co-(fac(2)-1)):fac(2)*co,(fac(3)*zz-(fac(3)-1)):fac(3)*zz, i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
    
elseif mod(fac(1),2) == 0 && mod(fac(2),2) == 0 && mod(fac(3),2) ~= 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):(fac(1)*ro+1),(fac(2)*co-(fac(2)-1)):(fac(2)*co+1),(fac(3)*zz-(fac(3)-1)):fac(3)*zz, i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
    
elseif mod(fac(1),2) == 0 && mod(fac(2),2) ~= 0 && mod(fac(3),2) == 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):(fac(1)*ro+1),(fac(2)*co-(fac(2)-1)):fac(2)*co,(fac(3)*zz-(fac(3)-1)):(fac(3)*zz+1), i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end

elseif mod(fac(1),2) ~= 0 && mod(fac(2),2) == 0 && mod(fac(3),2) == 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):fac(1)*ro,(fac(2)*co-(fac(2)-1)):(fac(2)*co+1),(fac(3)*zz-(fac(3)-1)):(fac(3)*zz+1), i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
    
elseif mod(fac(1),2) == 0 && mod(fac(2),2) == 0 && mod(fac(3),2) == 0
    for i = 1:S_csi(4)
        for ro = 1:S_csi(1)
            for co = 1:S_csi(2)
                for zz = 1:S_csi(3)
                    DP = csi_new((fac(1)*ro-(fac(1)-1)):(fac(1)*ro+1),(fac(2)*co-(fac(2)-1)):(fac(2)*co+1),(fac(3)*zz-(fac(3)-1)):(fac(3)*zz+1), i).*h1;
                    csi_sum(ro,co,zz,i) = sum(DP(:));
                end
            end
        end
    end
else
    disp('Input wrong upsample factor');
end

% h1 = fspecial3('gaussian', [5,5,3]);
% csi_new = zeros([120,120,21,512]);
% csi_new(:,:,1:20,:) = csi;
% csi_sum = zeros(S_csi);
% for i = 1:512
%     for ro = 1:S_csi(1)
%         for co = 1:S_csi(2)
%             for zz = 1:S_csi(3)
%                 DP = csi_new((5*ro-4):5*ro,(5*co-4):5*co,(2*zz-1):(2*zz+1),i).*h1;
%                 csi_sum(ro,co,zz,i) = sum(DP(:));
%             end
%         end
%     end
%     
% end

% %% PSF with Gaussian Kernel  (Need to adapt to more general case!)
% csi_new = zeros([93,93,90,512]);
% csi_new(1:92,1:92,1:90,:) = csi;
% csi_sum = zeros(S_csi);
% for i = 1:512
%     for ro = 1:S_csi(1)
%         for co = 1:S_csi(2)
%             for zz = 1:S_csi(3)
%                 DP = csi_new((2*ro-1):(2*ro+1),(2*co-1):(2*co+1),(5*zz-4):5*zz,i).*h1;
%                 csi_sum(ro,co,zz,i) = sum(DP(:));
%             end
%         end
%     end
%     
% end

% %% sum of subvoxels
% if(sum(fac) > numel(fac))
%     csi_sum = zeros(S_csi);
% 
%     for ro = 1:S_csi(1)
%         for co = 1:S_csi(2)
%             for zz = 1:S_csi(3)
%                 csi_sum(ro,co,zz,:) = sum(sum(sum(csi( (ro-1)*fac(1)+1:ro*fac(1), (co-1)*fac(2)+1:co*fac(2), (zz-1)*fac(3)+1:zz*fac(3),: ),1),2),3);
%             end
%         end
%     end
%     csi_sum = csi_sum / prod(fac);
% else
%     csi_sum = csi;
% end


%% for analyzing

% Inds_pos = zeros(64,64);
% Inds_neg = zeros(64,64);
% 
% for ro = 1:64
%     for co = 1:64
%         [~,Inds_pos(ro,co)] = max(real(fftshift(fft([squeeze(csi_sum(ro,co,:));zeros(8*2048,1)]))));
%         [~,Inds_neg(ro,co)] = max(real(fftshift(fft([squeeze(csi_sum_neg(ro,co,:));zeros(8*2048,1)]))));
% 
%     end
% end 
% 
% WP_pos = Inds_pos.*mask;
% WP_neg = Inds_neg.*mask;
% 
% WP_pos(WP_pos==0) = nan;
% WP_neg(WP_neg==0) = nan;
% 
% boxplot([reshape(WP_pos,1,[])',reshape(WP_neg,1,[])'])



% figure;
% subplot(121); imagesc(WP_pos./9)
% subplot(122); imagesc(WP_neg./9)
% toc;

%% cutting k_space back to original resolution
% kSp_csi_imp = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(csi_imp,1),2),[],1),[],2),1),2); % to kSpa
% 
% kSp_csi_imp_LD = kSp_csi_imp(end/2-S_csi(1)/2+1:end/2+S_csi(1)/2,end/2-S_csi(2)/2+1:end/2+S_csi(2)/2,:);
% csi_kspc =  fftshift(fftshift(fft(fft(ifftshift(ifftshift(kSp_csi_imp_LD,1),2),[],1),[],2),1),2); % to iSpace
% 






        
%%
% rr = 32;
% cc = 32;
% 
% figure;
% subplot(121); hold on
% for ro=1:fac 
%     for co =1:fac
%         plot(squeeze(real(fftshift(fft(csi_imp(rr+ro,cc+co,:),[],3),3))))
%     end
% end
% 
% subplot(122); hold on
% for ro=1:fac 
%     for co =1:fac
%         plot(squeeze(real(fftshift(fft(csi_zf(rr+ro,cc+co,:),[],3),3))))
%     end
% end
% %%
% figure; hold on
% plot(squeeze(sum(sum( real(fftshift(fft(csi_zf(rr+1:rr+ro,cc+1:cc+co,:),[],3),3)),1),2)));
% plot(squeeze(sum(sum( real(fftshift(fft(csi_imp(rr+1:rr+fac,cc+1:cc+co,:),[],3),3)),1),2)));

%% plotting subvoxel spectra
% 
% figure; hold on; for r = -5:6; for c = 1; plot(real(fftshift(fft(squeeze(csi_imp( (ro-1)*fac+r, (co-1)*fac+c, : )))))); pause ; end; end
% figure; hold on; for r = 1:3; for c = 1; plot(real(fftshift(fft(squeeze(csi_imp( (ro)*fac+r, (co)*fac+c, : )))))); pause ; end; end
% figure; hold on; for r = 1:3; for c = 1:3; plot(real(fftshift(fft(squeeze(csi_zf( (ro-1)*fac+r, (co-1)*fac+c, : )))))) ; end; end



% rows =S_csi(1);cols = S_csi(2);slices = S_csi(3);
% [X,Y,Z] = meshgrid(1:cols, 1:rows, 1:slices);
% [X2, Y2, Z2] = meshgrid(0.5:0.5:cols,0.5:0.5:rows, 0.2:0.2:slices);

% out = zeros(92,92,90,S_csi(4));%zeros(S_map);%
% for i = 1:S_csi(4)
% % out(:,:,:,i) = interp3(X, Y, Z, csi(:,:,:,i), X2, Y2, Z2, 'linear', 0);
% % out(:,:,:,i) = interp3(X, Y, Z, csi(:,:,:,i), X2, Y2, Z2, 'cubic', 0);
% out(:,:,:,i) = interp3(X, Y, Z, csi(:,:,:,i), X2, Y2, Z2, 'spline', 0);
% end
% csi = out;

% out = zeros(92,92,90,S_csi(4));
% Up_factor = [2,2,9];%[2,2,5];
% for k = 1:S_csi(4)
%   if abs(csi(:,:,:,k)) == 0
%       out(:,:,:,k) = 0;
% %       out(:,:,:,k) = interp3(X, Y, Z, csi(:,:,:,k), X2, Y2, Z2, 'spline', 0);
%   else
%       csi_real = NLMUpsample2(real(csi(:,:,:,k)),Up_factor);  
%       csi_imag = NLMUpsample2(imag(csi(:,:,:,k)),Up_factor);
%       out(:,:,:,k) = csi_real + 1i*csi_imag;
%   end
% end
% csi = out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% Interpolation in Image Domain
% Upfac = round(S_map./S_csi(1:3)); 
% ny=Upfac(1)*S_csi(1);nx=Upfac(2)*S_csi(2);nz=Upfac(3)*S_csi(3); %% desired output dimensions
% [y1, x1, z1] = ndgrid(linspace(1,size(csi,1),ny), linspace(1,size(csi,2),nx), linspace(1,size(csi,3),nz));
% out = zeros([size(B0map),S_csi(4)]);
% Im_LowRes = out;
 
% for i = 1:S_csi(4)
% % out(:,:,:,i) = interp3(csi(:,:,:,i), x1, y1, z1, 'linear', 0);
% % out(:,:,:,i) = interp3(csi(:,:,:,i), x1, y1, z1, 'cubic', 0);
% % out(:,:,:,i) = interp3(csi(:,:,:,i), x1, y1, z1, 'spline', 0);
% out(:,:,:,i) = interp3(csi(:,:,:,i), x1, y1, z1, 'nearest', 0);
% end
% csi = out;



% %% Nonlocal Upsampling
% for k = 1:S_csi(4)
%   if abs(csi(:,:,:,k)) == 0
%       out(:,:,:,k) = interp3(csi(:,:,:,k), x1,y1,z1,'spline', 0);
%   else
%       csi_real = NLMUpsample2(real(csi(:,:,:,k)),Upfac);  
%       csi_imag = NLMUpsample2(imag(csi(:,:,:,k)),Upfac);
%       out(:,:,:,k) = csi_real + 1i*csi_imag;
%   end
% end
% csi = out;



% %% Total variation based upsampling
% tveps = 1e-4;
% lambdaTV = 1e-5;
% niter = 50;
% h = 1;
% dt = 0.1;
% for k = 1:S_csi(4)
%   if abs(csi(:,:,:,k)) == 0
%       out(:,:,:,k) = interp3(csi(:,:,:,k), x1,y1,z1,'spline', 0);
%       Im_LowRes(:,:,:,k) = out(:,:,:,k);
%   else
%       g_real = InitialInterpolation(real(csi(:,:,:,k)),[3,3,4]);
%       g_imag = InitialInterpolation(imag(csi(:,:,:,k)),[3,3,4]);
%       out(:,:,:,k) = g_real + 1i*g_imag;
%       Im_LowRes(:,:,:,k) = out(:,:,:,k);
%       for j = 1:niter
%           aa = gauss3filter(Im_LowRes(:,:,:,k),h);
%           aa0_real = my_downsample(real(aa),[3,3,4]);
%           aa0_imag = my_downsample(imag(aa),[3,3,4]);
%           aa_real = InitialInterpolation(aa0_real,[3,3,4]);
%           aa_imag = InitialInterpolation(aa0_imag,[3,3,4]);
%           aa = aa_real + 1i*aa_imag;
%           
%           temp1 = gauss3filter(out(:,:,:,k)-aa,h);
%           
%           
%           temp2_real = zeros(size(Im_LowRes(:,:,:,k)));
%           temp2_imag = zeros(size(Im_LowRes(:,:,:,k)));
%           temp2 = zeros(size(Im_LowRes(:,:,:,k)));
%           for i = 1:size(Im_LowRes(:,:,:,k),3)
%               [ux,uy,uxx,uyy,uxy] = gradSecondOrder(real(Im_LowRes(:,:,i,k)));
%               temp2_real(:,:,i) =lambdaTV* (uxx.*(uy.^2+tveps) - 2*uxy.*ux.*uy + uyy.*(ux.^2+tveps) ) ./ ((ux.^2+uy.^2+tveps).^(1.5)) ;
%               [vx,vy,vxx,vyy,vxy] = gradSecondOrder(imag(Im_LowRes(:,:,i,k)));
%               temp2_imag(:,:,i) =lambdaTV* (vxx.*(vy.^2+tveps) - 2*vxy.*vx.*vy + vyy.*(vx.^2+tveps) ) ./ ((vx.^2+vy.^2+tveps).^(1.5)) ;
%               temp2(:,:,i) = temp2_real(:,:,i) + 1i*temp2_imag(:,:,i);
%           end
%           
%           %     temp3 = reconstructImage_parallel(Im_LowRes, weights, weightsInx, f1*NeighSize, f2*NeighSize, f3*NeighSize1, Nthreads);
%           
%           Im_LowRes(:,:,:,k) = Im_LowRes(:,:,:,k) + dt*(temp1 + temp2);
%           
% %           Im_LowRes(Im_LowRes < 0) = 0;
%           
%       end
%   end       
% end
% csi = Im_LowRes;


function [bima]=InitialInterpolation(nima1,lf)

s=size(nima1).*lf;
ori=((1+lf)/2);

% reconstruc using spline interpolation
[x,y,z] = ndgrid(ori(1):lf(1):1-ori(1)+s(1),ori(2):lf(2):1-ori(2)+s(2),ori(3):lf(3):1-ori(3)+s(3));
[xi,yi,zi] = ndgrid(1:s(1),1:s(2),1:s(3));
bima2 = interpn(x,y,z,nima1,xi,yi,zi,'spline'); 

% deal with extreme slices
for i=1:floor(lf(1)/2)
  bima2(i,:,:) = bima2(floor(lf(1)/2)+1,:,:);
end
for i=1:floor(lf(2)/2)
  bima2(:,i,:) = bima2(:,floor(lf(2)/2)+1,:);
end
for i=1:floor(lf(3)/2)
  bima2(:,:,i) = bima2(:,:,floor(lf(3)/2)+1);
end

for i=1:floor(lf(1)/2)
  bima2(s(1)-i+1,:,:) = bima2(s(1)-floor(lf(1)/2),:,:);
end
for i=1:floor(lf(2)/2)
  bima2(:,s(2)-i+1,:) = bima2(:,s(2)-floor(lf(2)/2),:);  
end
for i=1:floor(lf(3)/2)
  bima2(:,:,s(3)-i+1) = bima2(:,:,s(3)-floor(lf(3)/2));
end

% mean correction
for i=1:lf(1):s(1)
for j=1:lf(2):s(2)
for k=1:lf(3):s(3)  
    tmp=bima2(i:i+lf(1)-1,j:j+lf(2)-1,k:k+lf(3)-1);  
    off=nima1((i+lf(1)-1)/lf(1),(j+lf(2)-1)/lf(2),(k+lf(3)-1)/lf(3))-mean(tmp(:));
    bima(i:i+lf(1)-1,j:j+lf(2)-1,k:k+lf(3)-1)=bima2(i:i+lf(1)-1,j:j+lf(2)-1,k:k+lf(3)-1)+off;
end
end
end

