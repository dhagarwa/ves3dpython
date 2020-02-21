function val = integratePatch(p, f) %p is patch, f is the function. f should be compactly supported
%Integration using trapezoidal rule
%f is a scalar function
    fj = f.*(p.J); %multipying by jacobian
    fjh = fj*p.h_u*p.h_v;
    val = sum(fjh);



end