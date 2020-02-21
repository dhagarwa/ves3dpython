function [] = wrapper()
    
%     p1 = standardSpherePatch(16, 16, 1, 1);
%     p2 = standardSpherePatch(16, 16, 2, 1);
%     p3 = standardSpherePatch(16, 16, 3, 1);
%     p4 = standardSpherePatch(16, 16, 4, 1);
%     p5 = standardSpherePatch(16, 16, 5, 1);
%     p6 = standardSpherePatch(16, 16, 6, 1);
    
    m = 128;
    n = 128;
    R = 1;
    val = 0;
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       %size(patch.x)
       f = foo(patch.nodes);
       val = val + integratePatch(patch, f);
        
    end

   (val - 4*pi*R^2)/(4*pi*R^2)
end

