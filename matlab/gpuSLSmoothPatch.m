function val = gpuSLSmoothPatch(y, p, q)
    %p.delta = 0.6*(p.h_u)^0.9; 
    intgd = gpuSLSmoothIntegrand(gpuArray(y),  gpuArray(p.r), gpuArray(q), gpuArray(p.delta));
    quad = gpuintegratePatchVec(p, intgd);
    val = quad ;
end
