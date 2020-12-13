function [] = testpou()
%Testing the surface deriv function

      clc

    m = 63;
    n = m;
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
    m = 2*(m+1)-1;
    n = m;
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

    
    %using fftInterpolate code
    p1_pou_mat = ves2fft(p1.pou,p1.Nu,p1.Nv);
    p1_pou_mat_up = fftInterpolate(p1_pou_mat,[p2.Nu+1,p2.Nv+1]);
    p1_pou_up = fft2ves(p1_pou_mat_up,p2.Nu+1,p2.Nv+1);
    %figure;
    %stem3(reshape(p1_pou_up,[p2.Nu,p2.Nv]));    
    
    pou_up = upsample2(p1.pou,p1);
    %figure;
    %stem3(reshape(pou_up,[p2.Nu,p2.Nv]));
    max(abs(pou_up-p2.pou))
    max(abs(p1_pou_up-p2.pou))
    
    
    
    
   

end