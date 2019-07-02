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

% %% Interpolation in Image Domain
% Upfac = round(S_map./S_csi(1:3)); 
% ny=Upfac(1)*S_csi(1);nx=Upfac(2)*S_csi(2);nz=Upfac(3)*S_csi(3); %% desired output dimensions
% [y1, x1, z1] = ndgrid(linspace(1,size(csi,1),ny), linspace(1,size(csi,2),nx), linspace(1,size(csi,3),nz));
% out = zeros([size(B0map),S_csi(4)]);
%  
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



%% Frequency Align Data
csi = csi .* exp(myrepmat(2*pi*1i*B0map,size(csi)) .* myrepmat(t,size(csi)));



% h1 = fspecial3('gaussian', [3,3,3]);

% for i = 1:512
%     csi(:,:,:,i) = convn(csi(:,:,:,i), h1, 'same');
% end


%% Boxcar with Gaussian Averaging (Need to adapt to more general case!)
h1 = fspecial3('gaussian', [3,3,5]);
csi_new = zeros([138,138,41,512]);
csi_new(:,:,1:40,:) = csi;
csi_sum = zeros(S_csi);
for i = 1:512
    for ro = 1:S_csi(1)
        for co = 1:S_csi(2)
            for zz = 1:S_csi(3)
                DP = csi_new((3*ro-2):3*ro,(3*co-2):3*co,(4*zz-3):(4*zz+1),i).*h1;
                csi_sum(ro,co,zz,i) = sum(DP(:));
            end
        end
    end
    
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

