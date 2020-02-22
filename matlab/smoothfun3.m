function val = smoothfun3(x)
    
    a = erf(x) ;
    b = (8*x.^6 - 36*x.^4 + 6*x.^2 + 9);
    c = - (2/9)*x.*b.*exp(-1*x.^2) *(1/sqrt(pi));
    val = a + c;

end