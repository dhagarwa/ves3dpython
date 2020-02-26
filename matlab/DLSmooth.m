function val = DLSmooth(y, qx0, S) % y is  target point
    val = 0;
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       f = patch.q_dl;
       val = val + DLSmoothPatch(y, qx0, patch, f);
        
        
    end
    
    
    %Hard coded for on surface singularity correction
    factor = 0.5; % 1 for inside, 0.5 for on, 0 for outside
    val = val + factor*qx0;






end