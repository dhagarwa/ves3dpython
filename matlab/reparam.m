function [] = reparam(S)

    true_area = S.getArea()
    true_vol = S.getVolume()

    in = S.getInertia();
    
    %reparameterization
    tau = 0.00001;
    T = 0.02; nt = T/tau; minT = T/5;
    %plotSurface(S);
    
%    g0 = getPatchwiseLapSphere(S);
%     scatter3(y0(:,1), y0(:,2), y0(:,3), 10, g0(:, 1) ); colorbar;
%     axis([-1 1 -1 1 -1 1]); xlabel('XX'); ylabel('YY'); zlabel('ZZ');
%     figure;
    E = []; ti = []; grad = [];  delta = []; delta_min = []; areas = []; vols = []; moi = []; com = [];
    t_ = 0;
    all_grad = [];
    disp('Starting reparameterization');
    for ii=0:nt
        
        
        y0 = S.getPosition();
        
        [E0, density] = getParamEnergyBiLapSphere(S);
        n0 = S.getNormals();
        
        g0 = getPatchwiseBiLapSphere(S);
        %g0_x = S.blendSurfaceFunction(g0(:,1));
        %g0_y = S.blendSurfaceFunction(g0(:,2));
        %g0_z = S.blendSurfaceFunction(g0(:,3));
        %g0 = [g0_x g0_y g0_z];        
        pg0 = g0 - repmat(dot(n0,g0,2),[1,3]).*n0;
        u0 = -pg0;
        
        max_grad_norm = max(vecnorm(u0, 2, 2));
        grad_norm_l2 = getGradientNormSphere(S, u0);
        u0 = u0/max_grad_norm;
        cts = getInnerProductBiLapSphere(S, -u0);
        ts = min(tau, getInnerProductBiLapSphere(S, -u0));
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

            ti = [ti t_];
            areas = [areas S.getArea() ];
            vols = [vols S.getVolume()];
            in = S.getInertia();
            moi = [moi in(1,1)];
            com = [com S.getCenter()];
%             %figure;
%             %plotSurface(S);
%             f1 = figure;
%             ha=axes;
%             ha.FontSize = 16;
%             %ax = gca;
%             %ax.FontSize = 16;
%             scatter3(y0(:,1), y0(:,2), y0(:,3), 20, density, 'filled'); colorbar; caxis([0 3]);
%             axis([-1 1 -1 1 -1 1]);
%             xlabel('XX'); ylabel('YY'); zlabel('ZZ');
%             
%             t = title(['XZ plane point distribution, time = ', num2str(t_)]);
%             t.FontSize = 16;
%             view(ha,[0,0])
%             
%             %saveas(f1,['results/xz_symmdt_bilap',num2str(ii)],'jpg');
%             %#projection on the Y-Z plane
%             f2 = figure;
%             ha=axes;
%             %ax = gca;
%             %ax.FontSize = 16;
%             scatter3(y0(:,1), y0(:,2), y0(:,3), 20, density, 'filled'); colorbar;caxis([0 3]);
%             axis([-1 1 -1 1 -1 1]);
%             xlabel('XX'); ylabel('YY'); zlabel('ZZ');
%             
%             t =title(['YZ plane point distribution, time = ', num2str(t_)]);
%             t.FontSize = 16;
%             view(ha,[90,0])
%             %saveas(f2,['results/yz_symmdt_bilap',num2str(ii)],'jpg');
%             %#projection on the X-Y plane
%             f3 = figure;
%             %ax = gca;
%             %ax.FontSize = 16;
%             ha=axes;
%             scatter3(y0(:,1), y0(:,2), y0(:,3), 20, density, 'filled'); colorbar;caxis([0 3]);
%             axis([-1 1 -1 1 -1 1]);
%             xlabel('XX'); ylabel('YY'); zlabel('ZZ');
%             
%             t = title(['XY plane point distribution, time = ', num2str(t_)]);
%             t.FontSize = 16;
%             view(ha,[0,90])
%             %saveas(f3,['results/xy_symmdt_bilap',num2str(ii)],'jpg');
%             %disp('hello');
            timestep = ts;
        end   
        %blend u0
%         u0_x = S.blendSurfaceFunction(u0(:,1));
%         u0_y = S.blendSurfaceFunction(u0(:,2));
%         u0_z = S.blendSurfaceFunction(u0(:,3));
%         u0_blend = [u0_x u0_y u0_z];
        if (grad_norm_l2/grad(1)) < 10^-2 || abs(ts*grad_norm_l2) < 10^-7
            disp('GRAD DECAY COMPLETE');
            break
        end
        S.updateSurface(u0, ts);
        S.updateStale();
        %y1 = S.getPosition();
        t_ = t_ + ts;
%         if max(vecnorm(y1-y0,2, 2)) < 0.00001
%             break
%         end
        

    end
%     plotSurface(S); figure
%     scatter3(y0(:,1), y0(:,2), y0(:,3));
%     axis([-1 1 -1 1 -1 1])
%     xlabel('XX'); ylabel('YY'); zlabel('ZZ');
%     err = max(abs(vecnorm(y0,2,2)-1));

%       E
%       grad 
%       delta
%       delta_min
%       
%       areas
%       vols
%       moi
%       com      
%       ti
      

end

