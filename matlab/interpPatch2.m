function val = interpPatch2(f, x, p)
    %Function to interpolate f given as a column vector on a patch
    % x is two column  (u & v) query nodes given on the patch in (0,2*pi)
    % this function arranges f in matrix form as required interpfft2
    % function. Also need to ensure u, v with -inf values i.e., not in patch
    % are not counted (this is rn taken care of in fftinterp2 but thats not a good idea). 
    % p is given patch object
    x1 = 2*x(:, 1) - pi; % shift values to [-pi,pi] domain
    c =  2*pi/(pi + 2*p.eps_strip) * x(:, 2) - pi^2/(pi + 2*p.eps_strip); %shift values to [-pi,pi] domain
    x1 = [x1 c];
    u = x1(:, 1);
    v = x1(:, 2);
    f_up = upsample2(f, p); %upsample f by 2 times upsampling in each direction
    f_mat = ves2fft(f_up, p.uf*(p.Nu+1)-1, p.uf*(p.Nv+1)-1);%upsampled f in matrix form
    %pad one more zero column and row 
    f_mat = [f_mat zeros(p.uf*(p.Nu+1),1)];
    f_mat = [f_mat; zeros(1,p.uf*(p.Nu+1)+1)];
    f_mat = f_mat';
    m = p.Nu; n = p.Nv;
    [u_grid,v_grid]=ndgrid(2*pi*(0:2*m+2)/(2*m+2),  (2*pi)*(0:2*n+2)/(2*n+2));
    u_grid = u_grid - pi; v_grid = v_grid - pi; %move to [-pi,pi]
    F = griddedInterpolant(u_grid,v_grid,f_mat,'spline');
    val = F(u,v);
    val(isnan(val))=0;
    
    
    


end