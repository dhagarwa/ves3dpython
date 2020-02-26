function val = SLSmoothPatch(y, p, q )

    intgd = SLSmoothIntegrand(y,  p.r, q, 0.8*(p.h_u)^0.7);
    quad = integratePatchVec(p, intgd);
    val = quad ;
end