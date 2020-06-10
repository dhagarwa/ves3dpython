function f_app = fftinterp2(x, f)
    %funtion for inteprolation of function on 2d domain using ffts
    %f is matrix of function values on a uniform periodic square domain
    %[-pi, pi]^2, in the order as required in fftpde code. u in rows, v cols
    % x is the two column matrix points where value needs to be
    % approximated
    %NEED TO ACCELERATE IT, directly pass fft
    n = size(x, 1);
    nu = size(f, 2);
    nv = size(f, 1);
    f_app = zeros(n, 1);
    for ii=1:n
        u = x(ii, 1);
        v = x(ii, 2);
        if(u == -inf || v == -inf)
           f_app(ii) = 0; % NOTE u = -100 or v = -100 denotes not patch value, so this function is specific to the application not a general interpolation function anymore
            
        else
        f_app_col = zeros(1, nu);
        for jj=1:nu
            f_col = f(:, jj);
            f_app_col(jj) = fftinterp1(v, f_col);
            
        end
        f_app(ii) = fftinterp1(u, f_app_col');
        
        end
    end
    f_app = real(f_app);
end