function val = integratePatchDomain(p, f) %p is patch, f is the function. 
%Integration using trapezoidal rule
%f is a scalar function

    
    fjh = f*p.h_u*p.h_v;
    val = sum(fjh);



end