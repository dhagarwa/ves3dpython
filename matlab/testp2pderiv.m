function [] = testp2pderiv()
%Testing the surface deriv function

      clc

    m = 127;
    n = 127;
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
    %testlinks(S);
    p2 = patches(3);
    p1 = patches(2);
    
    f = sin(2*p2.u).*sin(2*p2.v).*sin(2*p2.u).*sin(2*p2.v);
    f_du = 4*sin(2*p2.u).*cos(2*p2.u).*sin(2*p2.v).*sin(2*p2.v);
    f_dv = 4*sin(2*p2.u).*sin(2*p2.u).*sin(2*p2.v).*cos(2*p2.v);

    jj = p2.numPatch;
    xq = p1.links(:, 2*jj-1:2*jj) %query points from links
    
    f_du_interp = 4*sin(2*xq(:,1)).*cos(2*xq(:,1)).*sin(2*xq(:,2)).*sin(2*xq(:,2));
    f_dv_interp = 4*sin(2*xq(:,1)).*sin(2*xq(:,1)).*sin(2*xq(:,2)).*cos(2*xq(:,2));
   
    for ii=1:size(f_du_interp, 1)
        
       if xq(ii, 1) == -Inf || xq(ii, 2) == -Inf
           
          f_du_interp(ii) = 0;
          f_dv_interp(ii) = 0;
       end
    end
    
    [f_du_app, f_dv_app] = S.p2p_deriv(f, p1, p2);
    abs([f_dv_app f_dv_interp]);
    %abs(f_du_interp)
    %abs(f_dv_app - f_dv_interp)
    max(abs(f_du_app - f_du_interp))
    max(abs(f_dv_app - f_dv_interp))
    
    
    
    
    

end