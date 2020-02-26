function y  = fooVec(x)
    
    %y = [zeros(size(x, 1), 1), -x(:, 3), x(:, 2)];
    %y = x;
    %y = ones(size(x, 1), 3);
    y = x.*x;

end