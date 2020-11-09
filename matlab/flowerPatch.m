function p = flowerPatch(m, n, numPatch, R)


p = Patch(m, n, R, numPatch, []);

    if numPatch==1 
        patch_r = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];

        end
        p.r = patch_r;
        p.rb = patch_r;
    elseif numPatch==2
        patch_r = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];

        end
        p.r = patch_r;
        p.rb = patch_r;
    elseif numPatch==5
         patch_r = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];

        end
        p.r = patch_r;
        p.rb = patch_r;
    elseif numPatch==6
        patch_r = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];

        end
        p.r = patch_r;
        p.rb = patch_r;
    elseif numPatch==3
        patch_r = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];

        end
        p.r = patch_r;
        p.rb = patch_r;
    elseif numPatch==4
        patch_r = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];

        end
        p.r = patch_r;
        p.rb = patch_r;
    end

    
p = p.update();    %updating x, y, z from r and pou
p = p.updateStale();
%p.J = abs(R^2*sin(p.u)); %Setting jacobian determinant

end


