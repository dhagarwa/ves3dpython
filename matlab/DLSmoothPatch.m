function val = DLSmoothPatch(y, qx0, p, q )

    intgd = DLSmoothIntegrand(y, qx0, p.r, q, sqrt(p.h_u), 1);
    quad = integratePatchVec(p, intgd);
    val = quad ;
end