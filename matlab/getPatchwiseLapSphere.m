function val = getPatchwiseLapSphere(S)
%calculate laplacian in u-v space. [x_uu y_uu z_uu]


        y0 = S.getPosition();

        patch = S.patches;
        x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
        y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
        z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];

        
        lap_x = S.getSphereLaplacian(x(:));
        lap_y = S.getSphereLaplacian(y(:));
        lap_z = S.getSphereLaplacian(z(:));
        

        
        val = -[lap_x lap_y lap_z];


end