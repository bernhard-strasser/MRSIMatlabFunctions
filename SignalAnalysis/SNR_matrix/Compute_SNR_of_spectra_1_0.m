function SNR_matrix = Compute_SNR_of_spectra_1_0(CSI_matrix,water_frequency,dwelltime,vecSize)


% IF THE PATH OF DAT FILE IS PROVIDED, READ THIS DAT FILE
if(ischar(CSI_matrix))

    dat_file = CSI_matrix;
    clear CSI_matrix
    CSI_matrix = read_csi_dat_1_3(dat_file);
    
end


% if the matrix is still of size (1,ROW,COL,SLC,vecSize) then reshape it or display error meassage if something is wrong 
if(numel(size(CSI_matrix)) > 5)    
    display([ char(10) 'ERROR: DIMENSION DISMATCH. Matrix must be either of size (ROW,COL,SLC,vecSize) or of size (1,ROW,COL,SLC,vecSize)' char(10) 'Did you forget to sum all channels of your data?' char(10) ])
    %quit force
    
elseif(numel(size(CSI_matrix)) == 5)
    
    if(size(CSI_matrix,1)>1)
        display([ char(10) 'ERROR: DIMENSION DISMATCH. Matrix must be either of size (ROW,COL,SLC,vecSize) or of size (1,ROW,COL,SLC,vecSize)' char(10) 'Did you forget to sum all channels of your data?' char(10) ])
        %quit force
    else    
        CSI_matrix = reshape(squeeze(CSI_matrix),[size(CSI_matrix,2) size(CSI_matrix,3) size(CSI_matrix,4) size(CSI_matrix,5)]);
    end
        
end


CSI_matrix = fftshift(fft(CSI_matrix,size(CSI_matrix,4),4),4);

chemshift_vector = compute_chemshift_vector_1_0(water_frequency,dwelltime,vecSize);
chemshift_vector_metabolites = chemshift_vector(chemshift_vector > 0.2 & chemshift_vector < 4.5);



size(chemshift_vector)
size(chemshift_vector_metabolites)
size(squeeze(real(CSI_matrix(round(size(CSI_matrix,1)/2),round(size(CSI_matrix,2)/2),round(size(CSI_matrix,3)/2),:))))


figure
plot(chemshift_vector,squeeze(real(CSI_matrix(round(size(CSI_matrix,1)/2),round(size(CSI_matrix,2)/2),round(size(CSI_matrix,3)/2),:)))); 
set(gca,'XDir','reverse');

figure
plot(chemshift_vector_metabolites,squeeze(real(CSI_matrix(round(size(CSI_matrix,1)/2),round(size(CSI_matrix,2)/2),round(size(CSI_matrix,3)/2),size(chemshift_vector_metabolites,1))))); 
set(gca,'XDir','reverse');


SNR_matrix = 0.5;
