function [] = checkFDM()
    
    m = 31; n = 1;
    ooa = 12;
    [D,coef] = getNonCompactFDmatrix(m, pi/(m+1),n, ooa);
    
    f = @(x) sin(x);
    x = (1:m)*pi/(m+1);
    Df = @(x) cos(x);
    error = max(abs(Df(x)' - D*f(x)'))
    
    
end