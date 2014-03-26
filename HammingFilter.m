function matrix_hammingfiltered = Hamming_filtering_1_2(matrix_in,do_hamming_filt_along_dimensions)
    

matrix_in_kspace = matrix_in;
for hamming_dim = do_hamming_filt_along_dimensions
    matrix_in_kspace = fftshift(ifft(ifftshift(matrix_in_kspace,hamming_dim),size(matrix_in,hamming_dim),hamming_dim),hamming_dim);
end


hamming_repmatted = ones(size(matrix_in));



for hamming_dim = do_hamming_filt_along_dimensions                      % Hamming filter in each dimension seperately
    
    
    RepmatToSizeOfMatrixIn = size(matrix_in);                           % The Hamming-filter must have the same size as the matrix_in
    RepmatToSizeOfMatrixIn(hamming_dim) = 1;                            % Do not repmat in that dimension, in which hamming filtering is performed.
    


    %calculate Hamming filter
    hamming_1D = (hamming(size(matrix_in,hamming_dim)));

    if(hamming_dim == 1)
        hamming_1D_LeadingOnes = hamming_1D;
    else
        
        reshape_to = horzcat(ones([1 hamming_dim-1]), numel(hamming_1D));
        hamming_1D_LeadingOnes = reshape(hamming_1D,reshape_to);
    end
    
    hamming_repmatted = repmat(hamming_1D_LeadingOnes, RepmatToSizeOfMatrixIn) .* hamming_repmatted;



end



matrix_hammingfiltered = matrix_in_kspace .* hamming_repmatted;



for hamming_dim = do_hamming_filt_along_dimensions
    matrix_hammingfiltered = fftshift(fft(ifftshift(matrix_hammingfiltered,hamming_dim),size(matrix_in,hamming_dim),hamming_dim),hamming_dim);
end