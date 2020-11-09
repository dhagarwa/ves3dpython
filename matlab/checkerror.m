function err = checkerror(f, p)
    %Function to check error in FFT representation 
    % f is periodic zero on boundary smooth function 
    % p is patch 
    
    
    m = p.Nu; n = p.Nv; 
    
    f_mat = ves2fft(f, m, n);
    
    f_mat_app = real(ifft2(fft2(f_mat)));
    
    err = max(max(abs(f_mat - f_mat_app)));
    


end