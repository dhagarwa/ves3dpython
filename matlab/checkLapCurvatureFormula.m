function [] = checkLapCurvatureFormula()
%checking \Delta_{\gamma} r = -2*H \hat{n}
      clc

    m = 127;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    patches = [];
    for i=1:6
       patch =  standardEllipsoidPatch(m, n, i, R, 1, 0.4, 1);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = ones(size(patch.r,1), 3);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, 1, patch, patch.q_dl);
       
        
    end

    S = Surface(patches, [1, 0, 0]);
    %DLvalue = DLSmooth(trg,1,S)
    
    %assign tensile force as q_sl
    sigma = [(patches(1).x).^2 (patches(2).x).^2 (patches(3).x).^2 (patches(4).x).^2 (patches(5).x).^2 (patches(6).x).^2]; % tensile scalar
    sigma = sigma(:);
    fs = S.getTensileForce(sigma);
    for ii=1:6
        
        S.patches(ii).q_sl = fs((ii-1)*m*n+1:ii*m*n, :);
        
    end
    

    x = S.getPosition();
    lap_x = S.getLaplacian(x(:,1));
    lap_y = S.getLaplacian(x(:,2));
    lap_z = S.getLaplacian(x(:,3));
    [H, K] = S.getCurvature();
    N = S.getNormals();
    lap1 = [lap_x lap_y lap_z];
    lap2 = -2*repmat(H, [1,3]).*N;
    
    err = max(vecnorm(lap1-lap2, 2, 2))
   
    fb = S.getBendingForce();
    maxval1 = max(-H)
    maxval2 = max(H.^2)
    maxval3 = max(H.^2-K)
    maxval4 = max(-H.*(H.^2-K))
    maxval5 = max(S.getLaplacian(-H))
    maxval6 = max(S.getLaplacian(-H) + 2*(-H).*(H.^2-K))
    maxval7 = max(vecnorm(fb, 2, 2))
    maxval8 = max(lap_x)
end

