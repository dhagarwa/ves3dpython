function E0 = getParamEnergy(S)
%calculate mesh energy function
%1/2 * int_uv ((x_u^2 + x_v^2) + (y_u^2 + y_v^2) + (z_u^2 + z_v^2) duv )

        y0 = S.getPosition();

        patch = S.patches;
        x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
        y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
        z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];

        [x_du_app, x_dv_app] = S.patchwiseDerivFDM(x);
        [y_du_app, y_dv_app] = S.patchwiseDerivFDM(y);
        [z_du_app, z_dv_app] = S.patchwiseDerivFDM(z);


        x_du = x_du_app(:);
        y_du = y_du_app(:);
        z_du = z_du_app(:);

        x_dv = x_dv_app(:);
        y_dv = y_dv_app(:);
        z_dv = z_dv_app(:);
        
        intgd = x_du.^2 + x_dv.^2 + y_du.^2 + y_dv.^2 + z_du.^2 + z_dv.^2;
        E0 = 0;
        
        for ii=1:6
            p = patch(ii);
            m = p.Nu; n = p.Nv;
            
            E0 = E0 + 0.5*integratePatchDomain(p, intgd((ii-1)*m*n+1:ii*m*n));
            
            
        end


end