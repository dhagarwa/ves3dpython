function val = SLSmoothIntegrand(y,  x, q, delta)
    %Smoothened Stokes SL kernel for on surface evaluation
    %y is target point, x is source points vector 
    %On surface only
    % q is surface density vector
    
       %assume y is the target calculate potential on all the source points 

        y = repmat(y, size(x, 1), 1);
        r = (y - x);
        t1 = dot(r, q, 2);
        abs_r = sqrt(sum(r.^2,2));
        abs_r_delta = abs_r/delta;
        s2 = smoothfun2(abs_r_delta);
        abs_r(abs_r == 0) = 1; %to avoid division by zero
        t1r = s2./(abs_r.^3);
        t1 = t1.*t1r;
        val1 = diag(t1)*r;  
        
        s1 = smoothfun1(abs_r_delta);
        t2 = s1./(abs_r);
        val2 = diag(t2)*q;
        
        val = 1/(8*pi) * (val1 + val2);
        
        
end