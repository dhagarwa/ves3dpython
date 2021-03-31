function val = integrateSphere(S, f)
    val = 0;
  
    for ii=1:length(S.patches)
       patch = S.patches(ii);
       f_patch = f((ii-1)*patch.numNodes+1:ii*patch.numNodes);
       val = val + integrateSpherePatch(patch, f_patch);
        
        
    end
    
    



end