function val = interpPatch(f, x, p)
    %Function to interpolate f given as a column vector on a patch
    % x is two column  (u & v) query nodes given on the patch
    % this function arranges f in matrix form as required interpfft2
    % function. Also need to ensure u, v with -1 values i.e., not in patch
    % are not counted (this is rn taken care of in fftinterp2 but thats not a good idea). 
    % p is given patch
    x1 = 2*x - pi; % shift values to [-pi,pi] domain
    u = x1(:, 1);
    v = x1(:, 2);
    f_mat = ves2fft(f, p.Nu, p.Nv);
    val = fftinterp2([u(:), v(:)], f_mat);
    
    
    


end