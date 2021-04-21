function [] = gpuEvolveSurface2()
    %membrane properties
    kb = 0;
    Es = 2;
    Ed = 20;
    vesprops = [0;Es;Ed]; %kb, Es, Ed
    
    %flow properties
    shear_rate = 1;
    
    %Initialize a surface
    
    

    m = 63
    n = m;
    R = 1;
    

    %trg = [0.997425 -0.0507148 -0.0507148];
    %norm_trg = sqrt(norm(trg))
    patches = [];
    numPatches = 6;
    for i=1:numPatches
       patch =  standardSpherePatch(m, n, i, R);
       patches = [patches patch];
        
    end
    patch.r;
    S = Surface(patches, vesprops);
    S.updateSurface(2*ones(size(S.getPosition(), 1), 3), 1);
    S.updateStale();
    ini_center = S.getCenter()
%     r_before_blend = S.getPosition();
%     S.blendSurface();
%     r_after_blend = S.getPosition();
%     error_blend = norm(r_before_blend - r_after_blend)/(m*n*6)
    t = 0;
    dt = 0.0025
    tot_xdisp = 0;
    numGPU = gpuDeviceCount;
    parpool(numGPU)
 
    iterstart = tic; 
    for nts=1:70
        %Using patches directly
%         for i=1:numPatches
%            patch = S.patches(i);
%            patch = patch.updateStale();
%            patch.q_sl = patch.shearForce(Es, Ed);
%            u_inf_patch = bgshearFlow(shear_rate, patch.r);
%            u_inf = [u_inf; u_inf_patch]; 
%         end
%         
%         u1 = u_inf + SLSurface(S);
        
        %Using surface class
        u_inf = bgPoiseuilleFlow(shear_rate, S.getPosition());
        S.updateStale();
        S.interfacialForce();
	slstart = tic;
        u = u_inf + gpuSLSurface(S);
        gpuSLtime = toc(slstart)
        S.updateSurface(u, dt);
        S.updateStale();
         %r_before_blend = S.getPosition();
         %S.blendSurface();
         %r_after_blend = S.getPosition();
         %error_blend = norm(r_before_blend - r_after_blend)/(m*n*6)
        before_center = S.getCenter();
        tot_xdisp = tot_xdisp + abs(S.recenter_x())
        nts
        S.updateStale();
        after_center = S.getCenter()
	if mod(nts, 5) == 0
		y0 = S.getPosition();
		fi = figure;
		scatter3(y0(:, 1), y0(:, 2), y0(:, 3));
		saveas(fi, ['jobresults/scatter_64_bigt_P_', num2str(nts)], 'jpg');
	end        
    end
    y0 = S.getPosition();
%    scatter3(y0(:,1), y0(:,2), y0(:,3));
    %axis([-2 2 -2 2 -2 2])
    
%     figure
%     plotPatch(S, 1);
%     figure
%     plotPatch(S, 2);
%     figure
%     plotPatch(S, 3);
%     figure
%     plotPatch(S, 4);
%     figure
%     plotPatch(S, 5);
%     figure
%     plotPatch(S, 6);
%     
     gpuiter = toc(iterstart)
end
