function [] = testderiv()
%Testing the surface deriv function

      clc

    m = 31;
    n = 31;
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
    p = patches(6);
    
    f = sin(4*p.u).*sin(8*p.v);
    f_du = 4*cos(4*p.u).*sin(8*p.v);
    f_dv = 8*sin(4*p.u).*cos(8*p.v);

    f_du_up = 4*cos(4*p.u_up).*sin(8*p.v_up);
    f_dv_up = 8*sin(4*p.u_up).*cos(8*p.v_up);
    
    [f_du_app, f_dv_app] = p.grad_patch(f);
    
    max(abs(f_du_app - f_du))
    max(abs(f_dv_app - f_dv))
    
    
    
    
   

end