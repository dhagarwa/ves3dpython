function [] = SLTest()
      clc

    m = 63;
    n = 63;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
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
 
    
    
 
    val = SLSmooth(trg,  S)
    true_val = [0.287151, 0.191636, 0.191636];
    SLerror = norm(val - true_val)/norm(true_val)
   
   
 
end

