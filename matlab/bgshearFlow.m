function u_inf = bgshearFlow(shear_rate, x)

    u_inf = zeros(size(x, 1), 3);
    u_inf(:, 1) = shear_rate*x(:, 2);


end