function p = flowerPatch(m, n, numPatch, R)


p = Patch(m, n, R, numPatch, []);

if numPatch == 0
syms th ph
%th = u, ph = v
rho =  1+ exp(-3*(1/4)*sqrt(105/(2*pi))*sin(th)^2*cos(th)*cos(2*ph));
x = rho*sin(th)*cos(ph); y = rho*sin(th)*sin(ph); z = rho*cos(th); 

X  = [x;y;z]; 
Xu = (diff(X,th)); 
Xv = (diff(X,ph)); 

Xuu = (diff(X,th,2)); 
Xuv = (diff(Xv,th)); 
Xvv = (diff(X,ph,2)); 
disp('....');

E = (Xu.'*Xu);
F = (Xu.'*Xv);
G = (Xv.'*Xv);
W = sqrt(E*G - F^2);
disp('......');

nor = cross(Xu,Xv)/W;
disp('........');

L = (Xuu.'*nor);
M = (Xuv.'*nor);
N = (Xvv.'*nor);
disp('..........');

H = -(((E*N - 2*F*M + G*L))/W^2/2)
%K = simplify(((L*N - M^2)/W^2));
disp('............');

else 
    syms th ph
    H = 0*th*ph
end
    

    if numPatch==1 
        patch_r = [];
        H_true = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];
            th = u0; ph = v0;
            H_curr = 0;
            H_true = [H_true;H_curr];

        end
        p.r = patch_r;
        p.rb = patch_r;
        p.H_true = H_true;
    elseif numPatch==2
        patch_r = [];
        H_true = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];
            th = u0; ph = v0;
            H_curr = 0;
            H_true = [H_true;H_curr];
        end
        p.r = patch_r;
        p.rb = patch_r;
        p.H_true = H_true;
    elseif numPatch==5
         patch_r = [];
         H_true = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];
            th = u0; ph = v0;
            H_curr = 0;
            H_true = [H_true;H_curr];
        end
        p.r = patch_r;
        p.rb = patch_r;
        p.H_true = H_true;
    elseif numPatch==6
        patch_r = [];
        H_true = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];
            th = u0; ph = v0;
            H_curr = 0;
            H_true = [H_true;H_curr];
        end
        p.r = patch_r;
        p.rb = patch_r;
        p.H_true = H_true;
    elseif numPatch==3
        patch_r = [];
        H_true = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];
            th = u0; ph = v0;
            H_curr = 0;
            H_true = [H_true;H_curr];
        end
        p.r = patch_r;
        p.rb = patch_r;
        p.H_true = H_true;
    elseif numPatch==4
        patch_r = [];
        H_true = [];
        for ii=1:size(p.u, 1)
            u = p.u(ii); v = p.v(ii); 
            r = sph2cartPatch_mod(u,v,numPatch);
            
            ref_node = p2pmap_modified(r,numPatch,1);
            u0 = ref_node(1); v0 = ref_node(2);
            rho = flower_rho(u0,v0);
            patch_r = [patch_r; rho*sin(u0)*cos(v0), rho*sin(u0).*sin(v0), rho*cos(u0) ];
            th = u0; ph = v0;
            H_curr = 0;
            H_true = [H_true;H_curr];
        end
        p.r = patch_r;
        p.rb = patch_r;
        p.H_true = H_true;
    end

    
p = p.update();    %updating x, y, z from r and pou
p = p.updateStale();
%p.J = abs(R^2*sin(p.u)); %Setting jacobian determinant

end


