function [] = SLTest2()
      clc

    m = 31;
    n = 31;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg))
    patches = [];
    true_sol = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       sph_n = 2; sph_m = 1;
       patch.q_sl = sphVec(sph_n, sph_m, patch, patch.r);
       true_sol = [true_sol;sph_n/((2*sph_n+1)*(2*sph_n+3))*patch.q_sl];
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches, [0, 0, 0]);
 
    
    
 
    val = SLSurface(S);
    true_val = [0.287151, 0.191636, 0.191636];
    %SLerror = norm(val - true_val)/norm(true_val)
    SLerror = norm(val - true_sol)/norm(true_sol)
   
 
end

