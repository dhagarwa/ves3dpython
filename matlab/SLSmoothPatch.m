function val = SLSmoothPatch(y, p, q)
    %p.delta = 0.6*(p.h_u)^0.9; 
    intgd = SLSmoothIntegrand(y,  p.r, q, p.delta);
    quad = integratePatchVec(p, intgd);
    val = quad ;
end