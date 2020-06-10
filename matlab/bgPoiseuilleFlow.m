function u_inf = bgPoiseuilleFlow(alpha, x)

    u_inf = zeros(size(x, 1), 3);
    u_inf(:, 1) = alpha*(25 - (x(:,2) - 2).^2 - (x(:,3) - 2).^2);


end