function val = gpuintegratePatch(p, f) %p is patch, f is the function. d
%Integration using trapezoidal rule
%f is a scalar function
    f0 = gpuArray(f).*gpuArray(p.pou); %multiplying by partition of unity
    fj = f0.*gpuArray(p.J); %multipying by jacobian
    
    fjh = fj*gpuArray(p.h_u*p.h_v);
    val = sum(fjh);



end
