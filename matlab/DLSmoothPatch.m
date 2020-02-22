function val = DLSmoothPatch(y, p, q )

    intgd = DLSmoothIntegrand(y, p.r, q, sqrt(p.h_u), 0);
    intgd
    quad = integratePatchVec(p, intgd);
    val = quad ;
end