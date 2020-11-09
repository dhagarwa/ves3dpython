function r = flower_function(u, v)
    Ylm = real(compute_Ylm(3,2,u',v'));
    Ylm = Ylm';
    rho = 1 + exp(-3*Ylm);
    r = [rho.*sin(u).*cos(v) rho.*sin(u).*sin(v) rho.*cos(u)];
    


end