function Vnm = sphVec(n, m, p, r)
% n, m are degree of spherical harmonics, p is patch object
% r is Nx3 cartesian coordinate matrix 
% U is Nx2 matrix containing (u, v) in a row for N points of eval
    Ynm = zeros(size(r, 1), 1);
    for ii = 1:size(r, 1)        
        [u, v, rho] = cart2sph(r(ii, 1), r(ii, 2), r(ii, 3));        
        Ynm(ii) = compute_Ylm(n, m, u, v);

    end
    
    Vnm = p.surf_grad(Ynm) - (n+1)*[Ynm, Ynm, Ynm].*r;

    

end