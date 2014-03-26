function matrix_hammingfiltered = Hamming_filtering_1_1(matrix_in,do_hamming_filt_along_dimensions)



do_not_hamming_filt_along_dimensions = setxor(do_hamming_filt_along_dimensions,1:numel(size(matrix_in)));  % setxor: computed the "complement" of two vectors, eg setxor([2 3], [1 2 3 4 5]) = [1 4 5])
size_mat = size(matrix_in);                                                 % size_mat(do_not...) gives the size of the matrix along the dimensions, along no hamming_filt is done prod() gives the product of all these numbers
repmat_number = prod(size_mat(do_not_hamming_filt_along_dimensions));       % example: do_not...= [1 4 5], size_mat(do_not...) = [32 4 1024] (from matrix size(matrix_in) = [32 64 64 4 1024]); prod(...) = 32*4*1024
    
%calculate Hamming filter
hamming_ROW = (hamming(size(matrix_in,do_hamming_filt_along_dimensions(1))));
hamming_COL = (hamming(size(matrix_in,do_hamming_filt_along_dimensions(2))));
hamming2D = hamming_ROW*hamming_COL';
hamming_repmatted = repmat(hamming2D, [1 1 repmat_number]);
hamming_repmatted = reshape(hamming_repmatted,size(matrix_in));


% convert back to k-space
matrix_hammingfiltered = matrix_in;
for dimension = do_hamming_filt_along_dimensions
    matrix_hammingfiltered = fftshift(ifft(ifftshift(matrix_hammingfiltered,dimension),size(matrix_in,dimension),dimension),dimension);
end


%apply Hamming filter
matrix_hammingfiltered = matrix_hammingfiltered .* hamming_repmatted;


% convert back to image space
for dimension = do_hamming_filt_along_dimensions
    matrix_hammingfiltered = fftshift(fft(ifftshift(matrix_hammingfiltered,dimension),size(matrix_in,dimension),dimension),dimension);
end

