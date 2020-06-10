function val = interpPatch(f, x, p)
    %Function to interpolate f given as a column vector on a patch
    % x is two column  (u & v) query nodes given on the patch in (0,2*pi)
    % this function arranges f in matrix form as required interpfft2
    % function. Also need to ensure u, v with -1 values i.e., not in patch
    % are not counted (this is rn taken care of in fftinterp2 but thats not a good idea). 
    % p is given patch object
    x1 = 2*x(:, 1) - pi; % shift values to [-pi,pi] domain
    c =  2*pi/(pi + 2*p.eps_strip) * x(:, 2) - pi^2/(pi + 2*p.eps_strip); %shift values to [-pi,pi] domain
    x1 = [x1 c];
    u = x1(:, 1);
    v = x1(:, 2);
    f_mat = ves2fft(f, p.Nu, p.Nv);
    val = fftinterp2([u(:), v(:)], f_mat);
    
    
    


end