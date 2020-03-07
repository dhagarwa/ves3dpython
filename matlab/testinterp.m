function [] = testinterp(n)
%Test for 2D interpolation
    hu = 2*pi/n;
    hv = 2*pi/n;
    [u,v] = meshgrid(-pi:hu:pi-hu, -pi:hv:pi-hv);
    
    f = ftest(u, v);
    
    
    nf = 2*n;
    
    hu = 2*pi/nf;
    hv = 2*pi/nf;
    [u,v] = meshgrid(-pi:hu:pi-hu, -pi:hv:pi-hv);
    
    f_exact = ftest(u, v);
    f_app = fftinterp2([u(:), v(:)], f);
    norm(f_exact(:) - f_app)

end

function val = ftest(x, y)
    %val = exp(sin(x).^2 + cos(exp(cos(x).^5)) ) + exp(sin(x).^2).*exp(sin(y));
    val = sin(x).^2 + cos(y);
end