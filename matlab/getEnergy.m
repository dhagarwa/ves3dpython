function E0 = getEnergy(S)
%calculate  energy function - incorrect

        y0 = S.getPosition();
        grad_y0_x = S.getSurfaceGradient(y0(:,1));
        grad_y0_y = S.getSurfaceGradient(y0(:,2));
        grad_y0_z = S.getSurfaceGradient(y0(:,3));
        
        
        intgd = sum(grad_y0_x.^2,2) + sum(grad_y0_y.^2,2) + sum(grad_y0_z.^2,2);
        
        E0 = 0;
        
        for ii=1:6
            patch = S.patches(ii);
            m = patch.Nu; n = patch.Nv;
            E0 = E0 + integratePatch(patch, intgd((ii-1)*m*n+1:ii*m*n));
            
            
        end


end