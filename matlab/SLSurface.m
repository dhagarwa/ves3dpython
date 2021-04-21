function val = SLSurface(S)
    val = [];
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       for jj=1:size(patch.r, 1)
          y = patch.r(jj, :);
          pot = SLSmooth(y, S);
          val = [val;pot];
       end
        
    end



end
