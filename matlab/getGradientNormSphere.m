function E0 = getGradientNormSphere(S, g)
%calculate L^2 norm of gradient pf parameterization energy


        intgd = sum(g.^2, 2);
        E0 = integrateSphere(S, intgd);
        
        
        

end