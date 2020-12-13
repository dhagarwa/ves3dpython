function val = interpPatch4(f, x, p)
    %Function to interpolate f given as a column vector on a patch
    % x is two column  (u & v) query nodes given on the patch in (0,pi)
    % this function arranges f in matrix form 
    % p is given patch object
    %f is not periodic smooth with zero boundary
    m = p.Nu; n = p.Nv;
    f_mat = reshape(f, [m,n]);
    
    u = x(:,1); v = x(:,2);

    [u_grid,v_grid]=ndgrid(pi*(1:m)/(m+1),  (pi)*(1:n)/(n+1));
    
    F = griddedInterpolant(u_grid,v_grid,f_mat,'spline');
    val = F(u,v);
    val(isnan(val))=0;
    %val(val==inf)=0;
    %val(val==-inf)=0;
    
    


end