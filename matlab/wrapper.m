function [] = wrapper()
    
%     p1 = standardSpherePatch(16, 16, 1, 1);
%     p2 = standardSpherePatch(16, 16, 2, 1);
%     p3 = standardSpherePatch(16, 16, 3, 1);
%     p4 = standardSpherePatch(16, 16, 4, 1);
%     p5 = standardSpherePatch(16, 16, 5, 1);
%     p6 = standardSpherePatch(16, 16, 6, 1);
    
    m = 64;
    n = 64;
    R = 1;
    val = 0;
    theta = pi/2;
    phi = pi/3;
    trg = 1*[R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       patches = [patches patch];
       patch.q_dl = fooVec(patch.r);
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches);
    val = DLSmooth(trg, S );
    %val = val + 0.5*fooVec([R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]);

   %(val - 4*pi*R^2)/(4*pi*R^2)
   error = norm(val - 0.5*[0, -trg(3), trg(2)]) 
end

