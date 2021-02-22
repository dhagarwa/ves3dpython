function [] = reparameterization()
      clc

    m = 15;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    patches = [];
    for i=1:6
       patch =  standardSphereSkewPatch(m, n, i, R, 4);
       
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
    
    area_err = abs(S.getArea() - 4*pi)/(4*pi)
    vol_err = abs(S.getVolume() - 4*pi/3)/(4*pi/3)
    patchint_err = abs(2*integratePatch_nopou(S.patches(1), ones(m*n, 1)) - 4*pi)/(4*pi)
    %reparameterization
    tau = 0.001;
    T = 0.02; nt = T/tau; minT = T/5;
    %plotSurface(S);
    y0 = S.getPosition();
    scatter3(y0(:,1), y0(:,2), y0(:,3));
    axis([-1 1 -1 1 -1 1])
    figure;
    E = []; t = []; grad = []; SLError = []; delta = []; areas = []; vols = []; moi = []; com = [];
    for ii=0:nt
        t_ = ii*tau;
        y0 = S.getPosition();
        iter = ii;
        E0 = getParamEnergy(S)

        n0 = S.getNormals();
        
        g0 = getPatchwiseLap(S);
        %g0_x = S.blendSurfaceFunction(g0(:,1));
        %g0_y = S.blendSurfaceFunction(g0(:,2));
        %g0_z = S.blendSurfaceFunction(g0(:,3));
        %g0 = [g0_x g0_y g0_z];        
        pg0 = g0 - repmat(dot(n0,g0,2),[1,3]).*n0;
        u0 = -pg0;
        
        max_grad_norm = max(vecnorm(u0, 2, 2))
        if mod(t_/minT,1) == 0
            E = [E E0];
            grad = [grad max_grad_norm];
            delta = [delta S.patches(1).delta];
            sle = SLTestreparam(S)
            SLError = [SLError sle];
            t = [t t_];
            areas = [areas abs(S.getArea() - 4*pi)/(4*pi)];
            vols = [vols abs(S.getVolume() - 4*pi/3)/(4*pi/3)];
            in = S.getInertia();
            moi = [moi abs(in(1,1)-8*pi/15)/(8*pi/15)];
            com = [com S.getCenter()];
        end   
        %blend u0
%         u0_x = S.blendSurfaceFunction(u0(:,1));
%         u0_y = S.blendSurfaceFunction(u0(:,2));
%         u0_z = S.blendSurfaceFunction(u0(:,3));
%         u0_blend = [u0_x u0_y u0_z];
        
        S.updateSurface(u0, tau);
        S.updateStale();
        y1 = S.getPosition();
%         if max(vecnorm(y1-y0,2, 2)) < 0.00001
%             break
%         end
        

    end
    %plotSurface(S);
    scatter3(y0(:,1), y0(:,2), y0(:,3));
    axis([-1 1 -1 1 -1 1])
    err = max(abs(vecnorm(y0,2,2)-1));
%     figure
%     plot(t, E, 'DisplayName', 'Energy');
%     title('Energy vs iterations')
%     figure
%     plot(t, grad, 'DisplayName', 'Gradient');
%     title('Gradient vs iterations')
%     figure 
%     plot(t, SLError, 'DisplayName', 'Single layer error');
%     title('Single layer error vs iterations')
%     figure 
%     plot( t, delta, 'DisplayName','Mesh Size');
%     title('Mesh size vs iterations')
      E
      grad 
      delta
      SLError
      areas
      vols
      moi
      com      
      t
      

end

