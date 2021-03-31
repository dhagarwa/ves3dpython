function [] = testnormalsFDM2Symm()
%Testing the surface deriv function

      clc

    m = 31;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg));
    patches = [];
    for i=1:6
       patch =  standardSpherePatchSymmetrical(m, n, i, R);
       
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
    
    %isnan(x_du_app)
    %isnan(x_dv_app)
    %isnan(y_du_app)
    %isnan(y_dv_app)
    %isnan(z_du_app)
    %isnan(z_dv_app)
    
    
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
    normals_fdm_new = D*temp_cross;
    
    nx = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
    ny = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
    nz = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];
    
    true_normals = [nx(:) ny(:) nz(:)];
    
    error_fdm_new = max(vecnorm(normals_fdm_new-true_normals, 2, 2))
    
    
    normals_fdm_new_blend_x = S.blendSurfaceFunction(normals_fdm_new(:,1));
    normals_fdm_new_blend_y = S.blendSurfaceFunction(normals_fdm_new(:,2));
    normals_fdm_new_blend_z = S.blendSurfaceFunction(normals_fdm_new(:,3));
    normals_fdm_new_blend = [normals_fdm_new_blend_x normals_fdm_new_blend_y normals_fdm_new_blend_z ];
    %isnan(normals_fdm_new)
    %calculate normals using fdm

    

    
    comp = [normals_fdm_new true_normals];
    
    comp_blend = [normals_fdm_new_blend true_normals];
    
    error_fdm_new_blend = max(vecnorm(normals_fdm_new_blend-true_normals, 2, 2))
    
    
      

end