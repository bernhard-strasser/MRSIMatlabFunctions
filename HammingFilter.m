function matrix_hammingfiltered = Hamming_filtering_1_0(matrix_in,filter_along_dimension)

% hamming_dimension = size(1,do_hamming_filt_along_dimensions);
% do_not_hamming_filt_along_dimensions = setxor(do_hamming_filt_along_dimensions,1:numel(size(matrix)));  % setxor: computed the "complement" of two vectors, eg setxor([1 2 5], [1 2 3 4 5 6 7 8]) = [3 4 6 7 8])
% 
% %hamming_nD = zeros(size(matrix));
% % convert back to k-space
% for dimension = do_hamming_filt_along_dimensions
%     matrix = fftshift(ifft(matrix,size(matrix,dimension),dimension));
% end
    


dim1 = filter_along_dimension(1);
dim2 = filter_along_dimension(2);
ROW = size(matrix_in,dim1);
COL = size(matrix_in,dim2);


%convert back to k-space
matrix_hammingfiltered = fftshift(fftshift(ifft(ifft(matrix_in,ROW,dim1),COL,dim2),dim1),dim2);

%calculate Hamming filter
hamming_ROW = (hamming(ROW));
hamming_COL = (hamming(COL));
hamming2D = hamming_ROW*hamming_COL';

%apply Hamming filter
matrix_hammingfiltered=matrix_hammingfiltered.*reshape(repmat(hamming2D,1,size(matrix_hammingfiltered,4)),ROW,COL,1,size(matrix_hammingfiltered,4));

%convert back to image space
matrix_hammingfiltered = fft2(ifftshift(ifftshift(matrix_hammingfiltered,1),2),ROW,COL);
    