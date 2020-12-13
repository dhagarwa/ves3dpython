function [] = testcurvatureFlowerFDM()
%Testing the surface curvature function for the sphere

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
       patch =  flowerPatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches, [0, 0, 0]);
    patch = patches;
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
    
    x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
    y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
    z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];
    
    
    [x_du_app, x_dv_app] = S.patchwiseDerivFDM(x);
    [y_du_app, y_dv_app] = S.patchwiseDerivFDM(y);
    [z_du_app, z_dv_app] = S.patchwiseDerivFDM(z);
    
    x_du = x_du_app(:);
    y_du = y_du_app(:);
    z_du = z_du_app(:);
    
    x_dv = x_dv_app(:);
    y_dv = y_dv_app(:);
    z_dv = z_dv_app(:);
    
    E = x_du.*x_du + y_du.*y_du + z_du.*z_du;
    F = x_du.*x_dv + y_du.*y_dv + z_du.*z_dv;
    G = x_dv.*x_dv + y_dv.*y_dv + z_dv.*z_dv;
    W = (E.*G - F.^2).^(0.5);
    
    
    r_du = [x_du y_du z_du];
    r_dv = [x_dv y_dv z_dv];
    
    temp_cross = cross(r_du, r_dv, 2);
    vecn = vecnorm(temp_cross, 2, 2).^(-1);
    D = spdiags(vecn(:),0,length(vecn),length(vecn));
    normals_fdm_new = D*temp_cross; %normals calculated
    %--------get true values on patch 1
    [H_true,n_true] = flowerTrueCurvature(m);
    n_true = reshape(n_true, [m*n, 3]);
    nx_true = n_true(:,1); ny_true = n_true(:,2); nz_true = n_true(:,3);
    %--------------
    normals_fdm_new_blend_x = S.blendSurfaceFunction(normals_fdm_new(:,1),nx_true);
    normals_fdm_new_blend_y = S.blendSurfaceFunction(normals_fdm_new(:,2),ny_true);
    normals_fdm_new_blend_z = S.blendSurfaceFunction(normals_fdm_new(:,3),nz_true);
    normals_fdm_new_blend = [normals_fdm_new_blend_x normals_fdm_new_blend_y normals_fdm_new_blend_z ];    
    
    [x_duu_app, x_duv_app] = S.patchwiseDerivFDM(x_du_app);
    [x_dvu_app, x_dvv_app] = S.patchwiseDerivFDM(x_dv_app);
    [y_duu_app, y_duv_app] = S.patchwiseDerivFDM(y_du_app);
    [y_dvu_app, y_dvv_app] = S.patchwiseDerivFDM(y_dv_app);    
    [z_duu_app, z_duv_app] = S.patchwiseDerivFDM(z_du_app);
    [z_dvu_app, z_dvv_app] = S.patchwiseDerivFDM(z_dv_app);
    
    size(x_du_app);
    size(x_duu_app);
    
    x_duu = x_duu_app(:);
    x_duv = x_duv_app(:);
    x_dvv = x_dvv_app(:);
    y_duu = y_duu_app(:);
    y_duv = y_duv_app(:);
    y_dvv = y_dvv_app(:);
    z_duu = z_duu_app(:);
    z_duv = z_duv_app(:);
    z_dvv = z_dvv_app(:);
    
    nx = normals_fdm_new_blend(:,1);
    ny = normals_fdm_new_blend(:,2);
    nz = normals_fdm_new_blend(:,3);
    
%     normals_fdm = [patch(1).nr;patch(2).nr;patch(3).nr;patch(4).nr;patch(5).nr;patch(6).nr];
%     nx = normals_fdm(:, 1);
%     ny = normals_fdm(:, 2);
%     nz = normals_fdm(:, 3);
    L = x_duu.*nx + y_duu.*ny + z_duu.*nz;
    M = x_duv.*nx + y_duv.*ny + z_duv.*nz;
    N = x_dvv.*nx + y_dvv.*ny + z_dvv.*nz;
    

    
    H = -0.5*(E.*N - 2*F.*M + G.*L)./ W.^2;
    H_blend = S.blendSurfaceFunction(H,H_true);
 
    
    error_blend = max(abs(H_blend(1:(m*n))-H_true))/max(abs(H_true))
    

    
    
   

end