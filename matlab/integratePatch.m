function val = integratePatch(p, f) %p is patch, f is the function. d
%Integration using trapezoidal rule
%f is a scalar function
    f0 = f.*(p.pou); %multiplying by partition of unity
    fj = f0.*(p.J); %multipying by jacobian
    
    fjh = fj*p.h_u*p.h_v;
    val = sum(fjh);



end