function[wdth]=op_CalcFWHMOfSpectra(InData,SearchPeakFromTo_ppm,Mask,DoSpectralFFT_flag)

[ppm_vector] = compute_chemshift_vector(InData);

if(~exist('DoSpectralFFT_flag','var'))
    DoSpectralFFT_flag = 0;
end
if(DoSpectralFFT_flag)
    InData.Data = fftshift(fft(InData.Data,[],4),4);
end
if(~exist('Mask','var') || isempty(Mask))
    Mask = ones(size_MultiDims(InData.Data,[1 2 3]));
end

%%
wdth = NaN(size_MultiDims(InData.Data,[1 2 3]));
for xVox = 1:size(InData.Data,1)
    for yVox = 1:size(InData.Data,2)
        for zVox = 1:size(InData.Data,3)

            if(Mask(xVox,yVox,zVox) == 0)
                wdth(xVox,yVox,zVox) = NaN;
                continue;
            end
            CurData = squeeze(InData.Data(xVox,yVox,zVox,:,:,:,:,:));
            if(sum(size(CurData) > 1) > 1)
                error('Error in op_CalcFWHMOfSpectra: Data is too high-dimensional. Rewrite function pleeease!')
            end
            if(~isreal(CurData))
                CurData = abs(CurData);
            end

            CurData_max=max(CurData);
            maxpos = find(CurData==CurData_max);
            % PeakRegion_LowerPt = find(min(abs(freq_vector - 1.9)) == abs(freq_vector - 1.9));
            PeakRegion_LowerPt = find(min(abs(ppm_vector - max(SearchPeakFromTo_ppm))) == abs(ppm_vector - max(SearchPeakFromTo_ppm)));
            PeakRegion_UpperPt = find(min(abs(ppm_vector - min(SearchPeakFromTo_ppm))) == abs(ppm_vector - min(SearchPeakFromTo_ppm)));


            try
                data1(:,1) = ppm_vector(PeakRegion_LowerPt:maxpos-1);
                data1(:,2) = CurData(PeakRegion_LowerPt:maxpos-1);
                data2(:,1) = ppm_vector(maxpos+1:PeakRegion_UpperPt);
                data2(:,2) = CurData(maxpos+1:PeakRegion_UpperPt);   
            catch wdth(xVox,yVox,zVox) = NaN;
            end    
        
        
            try
                x_halfmax1 = interp1(data1(:,2),data1(:,1),(CurData_max./2));
                x_halfmax2 = interp1(data2(:,2),data2(:,1),(CurData_max./2));
        
                wdth(xVox,yVox,zVox) = x_halfmax1-x_halfmax2;
            catch 
                wdth(xVox,yVox,zVox) = NaN;
            end


        end
    end
end
