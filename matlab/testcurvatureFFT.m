function [] = testcurvatureFFT()
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
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches, [0, 0, 0]);
    patch = patches;
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
    
    x = [patch(1).x patch(2).x patch(3).x patch(4).x patch(5).x patch(6).x];
    y = [patch(1).y patch(2).y patch(3).y patch(4).y patch(5).y patch(6).y];
    z = [patch(1).z patch(2).z patch(3).z patch(4).z patch(5).z patch(6).z];
    
    [x_du_app, x_dv_app] = S.derivFFT(x);
    [y_du_app, y_dv_app] = S.derivFFT(y);
    [z_du_app, z_dv_app] = S.derivFFT(z);
    
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
    normals_fft = D*temp_cross; %normals calculated
    
    
    
    [x_duu_app, x_duv_app] = S.derivFFT(x_du_app);
    [x_dvu_app, x_dvv_app] = S.derivFFT(x_dv_app);
    [y_duu_app, y_duv_app] = S.derivFFT(y_du_app);
    [y_dvu_app, y_dvv_app] = S.derivFFT(y_dv_app);    
    [z_duu_app, z_duv_app] = S.derivFFT(z_du_app);
    [z_dvu_app, z_dvv_app] = S.derivFFT(z_dv_app);
    
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
    
    nx = normals_fft(:,1);
    ny = normals_fft(:,2);
    nz = normals_fft(:,3);
    L = x_duu.*nx + y_duu.*ny + z_duu.*nz;
    M = x_duv.*nx + y_duv.*ny + z_duv.*nz;
    N = x_dvv.*nx + y_dvv.*ny + z_dvv.*nz;
    
    nx_t = [patch(1).x patch(2).x patch(3).x patch(4).x patch(5).x patch(6).x];
    ny_t = [patch(1).y patch(2).y patch(3).y patch(4).y patch(5).y patch(6).y];
    nz_t = [patch(1).z patch(2).z patch(3).z patch(4).z patch(5).z patch(6).z];
    
    true_normals = [nx_t(:) ny_t(:) nz_t(:)];
    
    error_fft = max(vecnorm(normals_fft-true_normals, 2, 2))
    
    H = 0.5*(E.*N - 2*F.*M + G.*L)./ W.^2;
    
    H_true = ones(size(H,1),1);

    
    error = max(abs(H-H_true))
    
    
    
    
    %max(abs(f_du_app - f_du))
    %max(abs(f_dv_app - f_dv))
    
    
    
    
   

end