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
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.902317 -0.39826 -0.164965];
    norm_trg = sqrt(norm(trg))
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       patch.q_sl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches);
    %val = DLSurface(S );
    %val = DLSmooth(trg, fooVec(trg), S);
    val = SLSmooth(trg,  S);
   val
%    tv = [];
%    for i=1:6
%        patch = S.patches(i);
%        true_val = 0.5*patch.q_dl ;
%        tv = [tv;true_val];
%    end
%    error = norm(val - tv)
end

