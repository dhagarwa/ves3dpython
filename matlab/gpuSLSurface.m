function val = gpuSLSurface(S)
    gpuval1 = zeros(length(S.patches)*size(S.patches(1).r, 1), 3, 'gpuArray');
  
    nodes = size(S.patches(1).r, 1);
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       gpuval2 = zeros(size(S.patches(1).r, 1), 3, 'gpuArray');
 
       parfor jj=1:size(patch.r, 1)
	  %gpuvaltemp = zeros(size(S.patches(1).r, 1), 3, 'gpuArray');
          y = patch.r(jj, :);
	  %gpuy = gpuArray(y);
          gpuval2(jj, :) = gpuSLSmooth(y, S);
          %val = [val;pot];
       end

       gpuval1((ii-1)*nodes+1:ii*nodes, :) = gpuval2;
   end
    val = gather(gpuval1);


end
