function [] = EvolveSurface()
    %membrane properties
    kb = 0;
    Es = 2;
    Ed = 20;
    vesprops = [0;Es;Ed]; %kb, Es, Ed
    
    %flow properties
    shear_rate = 5;
    
    %Initialize a surface
    
    clc

    m = 15;
    n = 15;
    R = 1;
    

    %trg = [0.997425 -0.0507148 -0.0507148];
    %norm_trg = sqrt(norm(trg))
    patches = [];
    numPatches = 6;
    for i=1:numPatches
       patch =  standardSpherePatch(m, n, i, R);
       patches = [patches patch];
        
    end
    S = Surface(patches, vesprops);
    S.updateSurface(2*ones(size(S.getPosition(), 1), 3), 1);
    S.updateStale();
    ini_center = S.getCenter()
%     r_before_blend = S.getPosition();
%     S.blendSurface();
%     r_after_blend = S.getPosition();
%     error_blend = norm(r_before_blend - r_after_blend)/(m*n*6)
    t = 0;
    dt = 0.0025;
    tot_xdisp = 0;
    tic 
    for nts=1:33
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
        u = u_inf + SLSurface(S);
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
        
    end
    
    figure
    plotPatch(S, 1);
    figure
    plotPatch(S, 2);
    figure
    plotPatch(S, 3);
    figure
    plotPatch(S, 4);
    figure
    plotPatch(S, 5);
    figure
    plotPatch(S, 6);
    
    toc
end