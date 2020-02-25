function val = DLSurface(S)
    val = [];
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       for jj=1:size(patch.r, 1)
          y = patch.r(jj, :);
          qx0 = patch.q_dl(jj, :);
          pot = DLSmooth(y, qx0, S);
          val = [val;pot];
       end
        
    end



end