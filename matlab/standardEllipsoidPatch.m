function p = standardEllipsoidPatch(m, n,  numPatch, R, a, b, c)
%ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
%a,b,c > 0

p = Patch(m, n, R, numPatch, []);

    if numPatch==1 
        p.r = [a*sin(p.u).*cos(p.v), b*sin(p.u).*sin(p.v), c*cos(p.u)];
        p.r_ = [a*sin(p.u_).*cos(p.v_), b*sin(p.u_).*sin(p.v_), c*cos(p.u_)];
        p.rb = [a*sin(p.U).*cos(p.V), b*sin(p.U).*sin(p.V), c*cos(p.U)];
    elseif numPatch==2
        %p.r = [R*sin(p.u).*cos(p.v), -R*sin(p.u).*sin(p.v), R*cos(p.u)];
        %before patch u-v direction alignment 
        
        p.r = [-a*sin(p.u).*cos(p.v), -b*sin(p.u).*sin(p.v), c*cos(p.u)];
        p.r_ = [-a*sin(p.u_).*cos(p.v_), -b*sin(p.u_).*sin(p.v_), c*cos(p.u_)];
        p.rb = [-a*sin(p.U).*cos(p.V), -b*sin(p.U).*sin(p.V), c*cos(p.U)];
    elseif numPatch==5
        p.r = [a*sin(p.u).*cos(p.v), -b*cos(p.u), c*sin(p.u).*sin(p.v)]; %u 0 to pi, v 0 to pi
        %p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)]; % u 0 to pi/2 , v 0 to 2*pi, patch is to be extended by mirror image
        p.rb = [a*sin(p.U).*cos(p.V), -b*cos(p.U), c*sin(p.U).*sin(p.V)];
        p.r_ = [a*sin(p.u_).*cos(p.v_), -b*cos(p.u_), c*sin(p.u_).*sin(p.v_)];
    elseif numPatch==6
        p.r = [a*sin(p.u).*cos(p.v), b*cos(p.u), -c*sin(p.u).*sin(p.v)]; %u 0 to pi, v 0 to pi
        %p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)]; % u pi/2 to pi, v 0 to 2*pi, patch is to be extended by mirror image
        p.rb = [a*sin(p.U).*cos(p.V), b*cos(p.U), -c*sin(p.U).*sin(p.V)];
        p.r_ = [a*sin(p.u_).*cos(p.v_), b*cos(p.u_), -c*sin(p.u_).*sin(p.v_)];
    elseif numPatch==3
        %p.r = [R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.r = [a*sin(p.u).*sin(p.v), -b*sin(p.u).*cos(p.v), c*cos(p.u)];
        p.rb = [a*sin(p.U).*sin(p.V), -b*sin(p.U).*cos(p.V), c*cos(p.U)];
        p.r_ = [a*sin(p.u_).*sin(p.v_), -b*sin(p.u_).*cos(p.v_), c*cos(p.u_)];
    elseif numPatch==4
        %p.r = [-R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.r = [-a*sin(p.u).*sin(p.v), b*sin(p.u).*cos(p.v), c*cos(p.u)];
        p.rb = [-a*sin(p.U).*sin(p.V), b*sin(p.U).*cos(p.V), c*cos(p.U)];
        p.r_ = [-a*sin(p.u_).*sin(p.v_), b*sin(p.u_).*cos(p.v_), c*cos(p.u_)];
    end

    
p = p.update();    %updating x, y, z from r and pou
p = p.updateStale();
%p.J = abs(R^2*sin(p.u)); %Setting jacobian determinant

end


