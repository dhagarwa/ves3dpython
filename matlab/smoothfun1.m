function val = smoothfun1(x)
    
    a = erf(x) ;
    b = (2*x.^2 - 5);
    c = - (2/3)*x.*b.*exp(-1*x.^2) *(1/sqrt(pi));
    val = a + c;

    %val = ones(size(x, 1), 1); %uncomment only for inside or outside surface
end