function val = SLSmooth(y, S) % y is target point
    val = 0;
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       f = patch.q_sl;
       val = val + SLSmoothPatch(y, patch, f);
        
        
    end
    
    



end