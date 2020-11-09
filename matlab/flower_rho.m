function rho = flower_rho(u, v)
    Ylm = real(compute_Ylm(3,2,u',v'));
    Ylm = Ylm';
    rho = 1 + exp(-3*Ylm);
    
    


end