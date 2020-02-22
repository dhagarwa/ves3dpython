function [] = wrapper()
    
%     p1 = standardSpherePatch(16, 16, 1, 1);
%     p2 = standardSpherePatch(16, 16, 2, 1);
%     p3 = standardSpherePatch(16, 16, 3, 1);
%     p4 = standardSpherePatch(16, 16, 4, 1);
%     p5 = standardSpherePatch(16, 16, 5, 1);
%     p6 = standardSpherePatch(16, 16, 6, 1);
    
    m = 16;
    n = 16;
    R = 1;
    val = 0;
    theta = 0.5;
    phi = 0.5;
    trg = [1+sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       %size(patch.x)
       f = fooVec(patch.nodes);
       %val = val + integratePatch(patch, f);
       val = val + DLSmoothPatch(trg, patch, f);
        
    end
    val = val + 0.5*[1 1 1];

   %(val - 4*pi*R^2)/(4*pi*R^2)
   val
end

