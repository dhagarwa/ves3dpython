function [] = testFDM()
    %Testing the surface deriv function

      clc

    m = 15;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg));
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches, [0, 0, 0]);
    p = patches(1);
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
    
    f = ones((m+2)*(n+2), 6);
    
    
    [f_du_app, f_dv_app] = S.derivFDM(f);
    
    (max(f_du_app))
    (max(f_dv_app))
    %max(abs(f_du_app - f_du))
    %max(abs(f_dv_app - f_dv))
    
    

end