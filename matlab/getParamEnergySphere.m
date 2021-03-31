function [E0, intgd] = getParamEnergySphere(S)
%calculate mesh energy function
%1/2 * int_uv ((x_u^2 + x_v^2) + (y_u^2 + y_v^2) + (z_u^2 + z_v^2) duv )

        y0 = S.getPosition();

        patch = S.patches;
        x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
        y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
        z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];

        grad_x = S.getSphereGradient(x(:));
        grad_y = S.getSphereGradient(y(:));
        grad_z = S.getSphereGradient(z(:));
        
        intgd = sum(grad_x.^2 + grad_y.^2 + grad_z.^2, 2);
        E0 = integrateSphere(S, intgd);
        
        
        

end