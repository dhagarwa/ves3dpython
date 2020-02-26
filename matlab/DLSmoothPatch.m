function val = DLSmoothPatch(y, qx0, p, q )

    intgd = DLSmoothIntegrand(y, qx0, p.r, q, p.n, 0.8*(p.h_u)^0.7, 1);
    quad = integratePatchVec(p, intgd);
    val = quad ;
end