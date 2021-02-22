function [] = testunique_uv_second_deriv()
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
       patch =  standardSphereSkewPatch(m, n, i, R, 1);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches, [0, 0, 0]);
    patch = patches;
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
    %point = patch(5).r(12,:)
    %patch(1).links_(:,10)
    %stem3(reshape(patch(1).links_(:,10), [patch(1).Nu+2,patch(1).Nv+2]));
    x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
    y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
    z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];
    
    [x_du_app, x_dv_app] = S.patchwiseDerivFDM(x);
    [y_du_app, y_dv_app] = S.patchwiseDerivFDM(y);
    [z_du_app, z_dv_app] = S.patchwiseDerivFDM(z);

    [x_duu_app, x_duv_app] = S.patchwiseDerivFDM(x_du_app);
    [x_dvu_app, x_dvv_app] = S.patchwiseDerivFDM(x_dv_app);
    [y_duu_app, y_duv_app] = S.patchwiseDerivFDM(y_du_app);
    [y_dvu_app, y_dvv_app] = S.patchwiseDerivFDM(y_dv_app);
    [z_duu_app, z_duv_app] = S.patchwiseDerivFDM(z_du_app);
    [z_dvu_app, z_dvv_app] = S.patchwiseDerivFDM(z_dv_app);
            
    
    [x_du_app_new, x_dv_app_new] = S.unique_uv_deriv(x);
    [y_du_app_new, y_dv_app_new] = S.unique_uv_deriv(y);
    [z_du_app_new, z_dv_app_new] = S.unique_uv_deriv(z);

    [x_duu_app_new, x_dvv_app_new] = S.unique_uv_second_deriv(x);
    [y_duu_app_new, y_dvv_app_new] = S.unique_uv_second_deriv(y);
    [z_duu_app_new, z_dvv_app_new] = S.unique_uv_second_deriv(z);
    
    
    x_du = x_du_app(:);
    y_du = y_du_app(:);
    z_du = z_du_app(:);
    
    x_dv = x_dv_app(:);
    y_dv = y_dv_app(:);
    z_dv = z_dv_app(:);
    
    r_du = [x_du y_du z_du];
    r_dv = [x_dv y_dv z_dv];

    x_du_new = x_du_app_new(:);
    y_du_new = y_du_app_new(:);
    z_du_new = z_du_app_new(:);
    
    x_dv_new = x_dv_app_new(:);
    y_dv_new = y_dv_app_new(:);
    z_dv_new = z_dv_app_new(:);
    
    r_du_new = [x_du_new y_du_new z_du_new];
    r_dv_new = [x_dv_new y_dv_new z_dv_new];
    
    %err_x = max(abs(x_du_new - x_du))
    err_unique_uv_deriv = max(vecnorm(r_du_new-r_du, 2, 2))
    err_unique_uv_deriv = max(vecnorm(r_dv_new-r_dv, 2, 2))
    

    err_sec_xuu = max(vecnorm(x_duu_app_new - x_duu_app, 2, 2))
    err_sec_xvv = max(vecnorm(x_dvv_app_new - x_dvv_app, 2, 2))
    err_sec_yuu = max(vecnorm(y_duu_app_new - y_duu_app, 2, 2))
    err_sec_yvv = max(vecnorm(y_dvv_app_new - y_dvv_app, 2, 2))
    err_sec_zuu = max(vecnorm(z_duu_app_new - z_duu_app, 2, 2))
    err_sec_zvv = max(vecnorm(z_dvv_app_new - z_dvv_app, 2, 2))
        
    
end