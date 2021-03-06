function [] = wrapper()
      clc
%     p1 = standardSpherePatch(16, 16, 1, 1);
%     p2 = standardSpherePatch(16, 16, 2, 1);
%     p3 = standardSpherePatch(16, 16, 3, 1);
%     p4 = standardSpherePatch(16, 16, 4, 1);
%     p5 = standardSpherePatch(16, 16, 5, 1);
%     p6 = standardSpherePatch(16, 16, 6, 1);
    
    m = 15;
    n = 15;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 0.0338093 0.0632528];
    norm_trg = sqrt(norm(trg))
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches);
 
    
    
    %val = DLSurface(S );
    val = DLSmooth(trg, fooVec(trg), S);
    %val = SLSmooth(trg,  S);
    true_val = [0.151684, 0.169363, 0.167644]
    error = norm(val - true_val)/norm(true_val)
   
   
    patch = patches(1);
    %test u-v derivatives
    f = sin(patch.U) + cos(patch.V) + (patch.U).^8; 
    [du, dv] = patch.grad_FDM(f);
    
    error_du = norm(du - cos(patch.U) - 8*(patch.U).^7)
    error_dv = norm(dv + sin(patch.V))
    %test tangents and normals
    error_normal = norm(patch.n - patch.rb/patch.R)/ (patch.numNodes)
    fs = patch.shearForce(1, 0.5)
    
    
    
    
%   [du, dv] = S.deriv(sin(2*pi/(pi + 2*patch.eps_strip)*([patch.nodes(:, 2) patch.nodes(:, 2) patch.nodes(:, 2) patch.nodes(:, 2)] + patch.eps_strip)) );
%    [du, dv] = S.deriv(ones(size(patch.r, 1), 6));
%    du
%    norm(du)
   %dv
%    tv = [];
%    for i=1:6
%        patch = S.patches(i);
%        true_val = 0.5*patch.q_dl ;
%        tv = [tv;true_val];
%    end
%    error = norm(val - tv)
end

