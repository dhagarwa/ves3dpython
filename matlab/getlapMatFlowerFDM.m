function [] = getlapMatFlowerFDM()
%get surface laplacian matrix 

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
    
    
    surf_lap_mat = zeros(m*n*6,m*n*6);
    %given function 
    phi = @(x,y,z) zeros(size(x,1),size(x,2));
    lap_phi = @(x,y,z) zeros(size(x,1),size(x,2)); %given laplacian
    
    %%%% f_du and f_dv calculation 
    f = phi(x,y,z); %Given function values 
    true_lap = lap_phi(x,y,z);
    
    
    for ii=1:m*n*6

    
        col = floor(ii/(m*n))+1;
        row = mod(ii,m*n)+1;
        f(row,col) = 1;

        [f_du_app, f_dv_app] = S.patchwiseDerivFDM(f);

        f_du = f_du_app(:);
        f_dv = f_dv_app(:);




        % calculate ((E*f_dv - F*f_du)/W) and take v derivative 
        t1 = (E.*f_dv - F.*f_du)./W;
        t1_app = reshape(t1, [size(x,1),size(x,2)]);

        [t1_du_app, t1_dv_app] = S.patchwiseDerivFDM(t1_app);
        t1_du = t1_du_app(:);
        t1_dv = t1_dv_app(:);


        % calculate ((G*f_du - F*f_dv)/W) and take u derivative

        t2 = (G.*f_du - F.*f_dv)./W;
        t2_app = reshape(t2, [size(x,1),size(x,2)]);

        [t2_du_app, t2_dv_app] = S.patchwiseDerivFDM(t2_app);
        t2_du = t2_du_app(:);
        t2_dv = t2_dv_app(:);

        %calculate surf_lap

        surf_lap = (t1_dv + t2_du)./W;

        surf_lap_blend = S.blendSurfaceFunction(surf_lap);
        true_lap = true_lap(:);


        surf_lap_mat(:,ii) = surf_lap_blend;
        f(row,col) = 0;
    end
    
    stem3(surf_lap_mat);
    S = svd(surf_lap_mat);
    plot(1:size(S,1),S);
    cond_number = cond(surf_lap_mat)

end

