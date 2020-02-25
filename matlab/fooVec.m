function y  = fooVec(x)
    
    y = [zeros(size(x, 1), 1), -x(:, 3), x(:, 2)];



end