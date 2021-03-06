function [] = reparameterization_sphere()
      clc

    m = 15;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    patches = [];
    a = 1; b = 0.1; c= 0.1;
    for i=1:6
       patch =  standardEllipsoidPatchSymmetrical(m, n, i, R, a, b, c);
       
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
    ecc = sqrt(1 - c^2/a^2); %eccentricity if c < a
    p = 1.6075;
    true_area = 4*pi* ((a^p*b^p + b^p*c^p + c^p*a^p)/3)^(1/p); true_vol = 4/3*pi*a*b*c; true_in_xx = 1/5*true_vol*(b^2 + c^2);
    area_err = abs(S.getArea() - true_area)/(true_area)
    vol_err = abs(S.getVolume() - true_vol)/(true_vol)
    in = S.getInertia();
    moi_err = abs(in(1, 1) - true_in_xx)/(true_in_xx)
    patchint_err = abs(2*integratePatch_nopou(S.patches(1), ones(m*n, 1)) - 4*pi)/(4*pi)
    %reparameterization
    tau = 0.0001;
    T = 1; nt = T/tau; minT = T/5;
    %plotSurface(S);
    y0 = S.getPosition(); 
%    g0 = getPatchwiseLapSphere(S);
%     scatter3(y0(:,1), y0(:,2), y0(:,3), 10, g0(:, 1) ); colorbar;
%     axis([-1 1 -1 1 -1 1]); xlabel('XX'); ylabel('YY'); zlabel('ZZ');
%     figure;
    E = []; ti = []; grad = []; SLError = []; delta = []; delta_min = []; areas = []; vols = []; moi = []; com = [];
    t_ = 0;
    all_grad = [];
    for ii=0:nt
        ii;
        
        y0 = S.getPosition();
        iter = ii;
        [E0, density] = getParamEnergyBiLapSphere(S);
        n0 = S.getNormals();
        
        g0 = getPatchwiseBiLapSphere(S);
        %g0_x = S.blendSurfaceFunction(g0(:,1));
        %g0_y = S.blendSurfaceFunction(g0(:,2));
        %g0_z = S.blendSurfaceFunction(g0(:,3));
        %g0 = [g0_x g0_y g0_z];        
        pg0 = g0 - repmat(dot(n0,g0,2),[1,3]).*n0;
        u0 = -pg0;
        
        max_grad_norm = max(vecnorm(u0, 2, 2))
        grad_norm_l2 = getGradientNormSphere(S, u0)
        u0 = u0/max_grad_norm;
        cts = getInnerProductBiLapSphere(S, -u0)
        ts = min(tau, getInnerProductBiLapSphere(S, -u0))
        all_grad = [all_grad grad_norm_l2];
%         if ii > 10 && (all_grad(ii-4) - all_grad(ii+2)) < 10^-5
%             %grad_diff = (all_grad(ii+1) - all_grad(ii+2))
%             %disp('hello')
%             tau = tau/10
%             ts = min(tau, getInnerProductSphere(S, -u0))
%             c_timestep = cts
%             ts = cts
%             disp('reduced dt_max')
%         end
           
        
        if mod(ii,500) == 0
            E = [E E0];
            grad = [grad grad_norm_l2];
            delta = [delta max([S.patches(1).delta, S.patches(2).delta, S.patches(3).delta, S.patches(4).delta, S.patches(5).delta, S.patches(6).delta])];
            delta_min = [delta_min min([S.patches(1).delta_min, S.patches(2).delta_min, S.patches(3).delta_min, S.patches(4).delta_min, S.patches(5).delta_min, S.patches(6).delta_min])];
            sle = SLTestreparam(S)
            SLError = [SLError sle];
            ti = [ti t_];
            areas = [areas abs(S.getArea() - true_area)/(true_area)];
            vols = [vols abs(S.getVolume() - true_vol)/(true_vol)];
            in = S.getInertia();
            moi = [moi abs(in(1,1)- true_in_xx)/(true_in_xx)];
            com = [com S.getCenter()];
            %figure;
            %plotSurface(S);
            f1 = figure;
            ha=axes;
            ha.FontSize = 16;
            %ax = gca;
            %ax.FontSize = 16;
            scatter3(y0(:,1), y0(:,2), y0(:,3), 20, density, 'filled'); colorbar; caxis([0 3]);
            axis([-1 1 -1 1 -1 1]);
            xlabel('XX'); ylabel('YY'); zlabel('ZZ');
            
            t = title(['XZ plane point distribution, time = ', num2str(t_)]);
            t.FontSize = 16;
            view(ha,[0,0])
            
            %saveas(f1,['results/xz_symmdt_bilap',num2str(ii)],'jpg');
            %#projection on the Y-Z plane
            f2 = figure;
            ha=axes;
            %ax = gca;
            %ax.FontSize = 16;
            scatter3(y0(:,1), y0(:,2), y0(:,3), 20, density, 'filled'); colorbar;caxis([0 3]);
            axis([-1 1 -1 1 -1 1]);
            xlabel('XX'); ylabel('YY'); zlabel('ZZ');
            
            t =title(['YZ plane point distribution, time = ', num2str(t_)]);
            t.FontSize = 16;
            view(ha,[90,0])
            %saveas(f2,['results/yz_symmdt_bilap',num2str(ii)],'jpg');
            %#projection on the X-Y plane
            f3 = figure;
            %ax = gca;
            %ax.FontSize = 16;
            ha=axes;
            scatter3(y0(:,1), y0(:,2), y0(:,3), 20, density, 'filled'); colorbar;caxis([0 3]);
            axis([-1 1 -1 1 -1 1]);
            xlabel('XX'); ylabel('YY'); zlabel('ZZ');
            
            t = title(['XY plane point distribution, time = ', num2str(t_)]);
            t.FontSize = 16;
            view(ha,[0,90])
            %saveas(f3,['results/xy_symmdt_bilap',num2str(ii)],'jpg');
            %disp('hello');
            timestep = ts
        end   
        %blend u0
%         u0_x = S.blendSurfaceFunction(u0(:,1));
%         u0_y = S.blendSurfaceFunction(u0(:,2));
%         u0_z = S.blendSurfaceFunction(u0(:,3));
%         u0_blend = [u0_x u0_y u0_z];
        if (grad_norm_l2/grad(1)) < 10^-2 %|| ts < 10^-8
            disp('GRAD DECAY COMPLETE');
            break
        end
        S.updateSurface(u0, ts);
        S.updateStale();
        y1 = S.getPosition();
        t_ = t_ + ts;
%         if max(vecnorm(y1-y0,2, 2)) < 0.00001
%             break
%         end
        

    end
    plotSurface(S); figure
    scatter3(y0(:,1), y0(:,2), y0(:,3));
    axis([-1 1 -1 1 -1 1])
    xlabel('XX'); ylabel('YY'); zlabel('ZZ');
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
      delta_min
      SLError
      areas
      vols
      moi
      com      
      ti
      

end

