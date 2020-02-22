function val = DLSmoothIntegrand(y, x, q, delta, default)
    %Smoothened Stokes kernel for on surface evaluation
    %y is target point, x is source points vector 
    %On surface only
    % q is surface density vector
    
    if default == 1
       %assume y = x, calculate potential on all the source points 
        
        
    elseif default == 0   % y is given explicitly
        q_y = [1 1 1]; %for now to check scheme we take 1, but need to modify for automatic value calculation
        n = x; %for now we know the normal is x for unit sphere but need to write code for it.
        y = repmat(y, size(x, 1), 1);
        q_y = repmat(q_y, size(x, 1), 1);
        r = (y - x);
        q_tilde = q - q_y;
        t1 = dot(r, q_tilde, 2);
        t2 = dot(r, n, 2);
        abs_r = sqrt(sum(r.^2,2));
        abs_r_delta = abs_r/delta;
        s3 = smoothfun3(abs_r_delta);
        t4 = s3./(abs_r.^5);
        t12 = t1.*t2;
        t124 = t12.*t4;
        val = (-3/(4*pi))*diag(t124)*r;
        
        
    end

end