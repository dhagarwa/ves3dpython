function [] = testp2pderivFFT()
%Testing the surface deriv function

      clc

    m = 63;
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
    p = patches(5);
    p2 = 1;
    
    d1 = p.deriv_links(:,4*p2-3);
    d2 = p.deriv_links(:,4*p2-2);
    d3 = p.deriv_links(:,4*p2-1);
    d4 = p.deriv_links(:,4*p2-0);
    
    figure;
    stem3(reshape(d1,[p.Nu,p.Nv]));
    figure;
    stem3(reshape(d2,[p.Nu,p.Nv]));
    figure;
    stem3(reshape(d3,[p.Nu,p.Nv]));
    figure;
    stem3(reshape(d4,[p.Nu,p.Nv]));
    
    
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
    
    %f = ones(p.numNodes, 6);
    
    
    %[f_du_app, f_dv_app] = S.derivFFT(f);
    
    %max(max(f_du_app))
    %max(max(f_dv_app))
    %max(abs(f_du_app - f_du))
    %max(abs(f_dv_app - f_dv))
    
    
    
    
   

end