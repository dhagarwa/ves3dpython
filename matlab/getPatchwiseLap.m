function val = getPatchwiseLap(S)
%calculate laplacian in u-v space. [x_uu y_uu z_uu]


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

        [x_duu_app, x_duv_app] = S.patchwiseDerivFDM(x_du_app);
        [x_dvu_app, x_dvv_app] = S.patchwiseDerivFDM(x_dv_app);
        [y_duu_app, y_duv_app] = S.patchwiseDerivFDM(y_du_app);
        [y_dvu_app, y_dvv_app] = S.patchwiseDerivFDM(y_dv_app);
        [z_duu_app, z_duv_app] = S.patchwiseDerivFDM(z_du_app);
        [z_dvu_app, z_dvv_app] = S.patchwiseDerivFDM(z_dv_app);
        
%         [x_duu_app, x_dvv_app] = S.unique_uv_second_deriv(x);
%         [y_duu_app, y_dvv_app] = S.unique_uv_second_deriv(y);
%         [z_duu_app, z_dvv_app] = S.unique_uv_second_deriv(z);
%            
        
        x_duu = x_duu_app(:);
        x_duv = x_duv_app(:);
        x_dvv = x_dvv_app(:);
        y_duu = y_duu_app(:);
        y_duv = y_duv_app(:);
        y_dvv = y_dvv_app(:);
        z_duu = z_duu_app(:);
        z_duv = z_duv_app(:);
        z_dvv = z_dvv_app(:);
        
        val = -[x_duu+x_dvv y_duu+y_dvv z_duu+z_dvv];


end