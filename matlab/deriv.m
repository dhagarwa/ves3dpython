function [f_du, f_dv] = deriv(f)
%derivative using ffts of function 2D function f on [-pi, pi]^2. 
% Nx is number of points on x-grid, Ny on y grid
% f should be Nx*Ny size matrix. 

    Nx = size(f, 1);
    Ny = size(f, 2);
    o = calc2([Nx, Ny]');
    
    o.use_hou_filtering;
    %o.use_twothirds_filtering;
    
    
    val =  o.grad(f);
    f_du = val(:, :, 1);
    f_dv = val(:, :, 2);
    
    


end