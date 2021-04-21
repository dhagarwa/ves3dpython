function psi = inclination(J)
%function to caluclate inclination from moment of inertia tensor


        [V,D] = eigs(J, 1, 'SM');
        if V(1) < 0 
            V = -V;
        end
        %V = getMin(V, D);
        psi = atan2(norm(cross(V,[1;0;0])),dot(V,[1;0;0]));
    




end
