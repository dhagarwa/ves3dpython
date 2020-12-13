function [] = testnormalsFFT()
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
    
    r_du = [x_du y_du z_du];
    r_dv = [x_dv y_dv z_dv];
    
    temp_cross = cross(r_du, r_dv, 2);
    vecn = vecnorm(temp_cross, 2, 2).^(-1);
    D = spdiags(vecn(:),0,length(vecn),length(vecn));
    normals_fft = D*temp_cross;
    
    %calculate normals using fdm
    normals_fdm = [patch(1).nr;patch(2).nr;patch(3).nr;patch(4).nr;patch(5).nr;patch(6).nr];
    
    nx = [patch(1).x patch(2).x patch(3).x patch(4).x patch(5).x patch(6).x];
    ny = [patch(1).y patch(2).y patch(3).y patch(4).y patch(5).y patch(6).y];
    nz = [patch(1).z patch(2).z patch(3).z patch(4).z patch(5).z patch(6).z];
    
    true_normals = [nx(:) ny(:) nz(:)];
    
    error_fft = max(vecnorm(normals_fft-true_normals, 2, 2))
    
    error_fdm = max(vecnorm(normals_fdm-true_normals, 2, 2))
    %max(abs(f_du_app - f_du))
    %max(abs(f_dv_app - f_dv))
    
    
    
    
   

end