function [] = wrapper()
    
%     p1 = standardSpherePatch(16, 16, 1, 1);
%     p2 = standardSpherePatch(16, 16, 2, 1);
%     p3 = standardSpherePatch(16, 16, 3, 1);
%     p4 = standardSpherePatch(16, 16, 4, 1);
%     p5 = standardSpherePatch(16, 16, 5, 1);
%     p6 = standardSpherePatch(16, 16, 6, 1);
    
    m = 8;
    n = 8;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = 1*[R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches);
    val = DLSurface(S );
    %val = val + 0.5*fooVec([R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]);

   %(val - 4*pi*R^2)/(4*pi*R^2)
   %error = norm(val - 0.5*[0, -trg(3), trg(2)]) 
   val
   tv = [];
   for i=1:6
       patch = S.patches(i);
       true_val = 0.5*patch.q_dl ;
       tv = [tv;true_val];
   end
   error = norm(val - tv)
end

