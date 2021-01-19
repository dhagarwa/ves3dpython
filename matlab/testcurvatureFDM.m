function [] = testcurvatureFDM()
%Testing the surface curvature function for the sphere

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
    
    normals_fdm_new_blend_x = S.blendSurfaceFunction(normals_fdm_new(:,1));
    normals_fdm_new_blend_y = S.blendSurfaceFunction(normals_fdm_new(:,2));
    normals_fdm_new_blend_z = S.blendSurfaceFunction(normals_fdm_new(:,3));
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
    
    nx_t = [patch(1).x patch(2).x patch(3).x patch(4).x patch(5).x patch(6).x];
    ny_t = [patch(1).y patch(2).y patch(3).y patch(4).y patch(5).y patch(6).y];
    nz_t = [patch(1).z patch(2).z patch(3).z patch(4).z patch(5).z patch(6).z];
    
    true_normals = [nx_t(:) ny_t(:) nz_t(:)];
    
    error_normals = max(vecnorm(normals_fdm_new_blend-true_normals, 2, 2))
    
    H = -0.5*(E.*N - 2*F.*M + G.*L)./ W.^2;
    H_blend = S.blendSurfaceFunction(H);
    
    H_true = ones(size(H,1),1);
    K_true = ones(size(H,1),1);
    
    error = max(abs(H-H_true))
    error_blend = max(abs(H_blend-H_true))
    %using surface class methods 
    n_blend2 = S.getNormals();
    [H_blend2, K_blend2] = S.getCurvature();
    fb = S.getBendingForce();
    error_nblend = max(vecnorm(n_blend2-true_normals, 2, 2))
    error_Hblend = max(abs(H_blend2-H_true))
    error_Kblend = max(abs(K_blend2-K_true))
    
    
    
%     syms u v
% 
% 
% x = sin(u)*cos(v); y = sin(u)*sin(v); z = cos(u); 
% 
% X  = [x;y;z]; 
% Xu = (diff(X,u)); 
% Xv = (diff(X,v)); 
% disp('..');
% 
% Xuu = (diff(X,u,2)); 
% Xuv = (diff(Xv,u)); 
% Xvv = (diff(X,v,2)); 
% disp('....');
% 
% E = (Xu.'*Xu);
% F = (Xu.'*Xv);
% G = (Xv.'*Xv);
% W = sqrt(E*G - F^2);
% disp('......');
% 
% nor = cross(Xu,Xv)/W;
% disp('........');
% 
% L = (Xuu.'*nor);
% M = (Xuv.'*nor);
% N = (Xvv.'*nor);
% disp('..........');
% 
% H = simplify(-((E*N - 2*F*M + G*L))/W^2/2)
% %K = ((L*N - M^2)/W^2);
% disp('............');
% u = [pi/2,pi/4]; v = [pi/2,pi/4];
% subs(H)
    
    %max(abs(f_du_app - f_du))
    %max(abs(f_dv_app - f_dv))
    
    
    
    
   

end