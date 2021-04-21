function val = gpuSLSmooth(y, S) % y is target point
    val = 0;
    gpuval = zeros(length(S.patches), 3, 'gpuArray');
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       f = patch.q_sl;
       gpuval(ii, :) = gpuSLSmoothPatch(y, patch, f);
        
        
    end
    for ii=1:length(S.patches)

	val = val + gpuval(ii, :);
    end
   
    



end
