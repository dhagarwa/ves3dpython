function val = SLSmoothPatch(y, p, q )

    intgd = SLSmoothIntegrand(y,  p.r, q, 5*(p.h_u)^0.9);
    quad = integratePatchVec(p, intgd);
    val = quad ;
end