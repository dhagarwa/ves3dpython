function [] = test_deriv()
    clc;
    Nx = 16;
    Ny = Nx;
    o = calc2([Nx, Ny]');
    
    o.use_hou_filtering;
    %o.use_twothirds_filtering;
    X=o.regulargrid
    
    f = sin(X(:, :, 1))
    size(f)
    %true_grad = ;
    approx_grad = o.grad(f)
    %error = norm(o.C(true_grad-approx_grad))


end