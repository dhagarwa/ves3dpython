function [] = testpou()
%Testing the surface deriv function

      clc

    m = 63;
    n = 63;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg));
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S1 = Surface(patches, [0, 0, 0]);
    p1 = patches(1);
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
 
    stem3(reshape(p1.pou,[p1.Nu,p1.Nv]));
    m = 127;
    n = 127;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg));
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S2 = Surface(patches, [0, 0, 0]);
    p2 = patches(1);

    
    
    
    pou_up = upsample2(p1.pou,p1);
    figure;
    stem3(reshape(pou_up,[p2.Nu,p2.Nv]));
    max(abs(pou_up-p2.pou))
    
    
    
    
   

end