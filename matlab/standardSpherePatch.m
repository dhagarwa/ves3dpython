function p = standardSpherePatch(m, n, numPatch, R)


p = Patch(m, n, R, numPatch, []);

    if numPatch==1 
        p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)];
        p.r_ = [R*sin(p.u_).*cos(p.v_), R*sin(p.u_).*sin(p.v_), R*cos(p.u_)];
        p.rb = [R*sin(p.U).*cos(p.V), R*sin(p.U).*sin(p.V), R*cos(p.U)];
    elseif numPatch==2
        %p.r = [R*sin(p.u).*cos(p.v), -R*sin(p.u).*sin(p.v), R*cos(p.u)];
        %before patch u-v direction alignment 
        
        p.r = [-R*sin(p.u).*cos(p.v), -R*sin(p.u).*sin(p.v), R*cos(p.u)];
        p.r_ = [-R*sin(p.u_).*cos(p.v_), -R*sin(p.u_).*sin(p.v_), R*cos(p.u_)];
        p.rb = [-R*sin(p.U).*cos(p.V), -R*sin(p.U).*sin(p.V), R*cos(p.U)];
    elseif numPatch==5
        p.r = [R*sin(p.u).*cos(p.v), -R*cos(p.u), R*sin(p.u).*sin(p.v)]; %u 0 to pi, v 0 to pi
        %p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)]; % u 0 to pi/2 , v 0 to 2*pi, patch is to be extended by mirror image
        p.rb = [R*sin(p.U).*cos(p.V), -R*cos(p.U), R*sin(p.U).*sin(p.V)];
        p.r_ = [R*sin(p.u_).*cos(p.v_), -R*cos(p.u_), R*sin(p.u_).*sin(p.v_)];
    elseif numPatch==6
        p.r = [R*sin(p.u).*cos(p.v), R*cos(p.u), -R*sin(p.u).*sin(p.v)]; %u 0 to pi, v 0 to pi
        %p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)]; % u pi/2 to pi, v 0 to 2*pi, patch is to be extended by mirror image
        p.rb = [R*sin(p.U).*cos(p.V), R*cos(p.U), -R*sin(p.U).*sin(p.V)];
        p.r_ = [R*sin(p.u_).*cos(p.v_), R*cos(p.u_), -R*sin(p.u_).*sin(p.v_)];
    elseif numPatch==3
        %p.r = [R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.r = [R*sin(p.u).*sin(p.v), -R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.rb = [R*sin(p.U).*sin(p.V), -R*sin(p.U).*cos(p.V), R*cos(p.U)];
        p.r_ = [R*sin(p.u_).*sin(p.v_), -R*sin(p.u_).*cos(p.v_), R*cos(p.u_)];
    elseif numPatch==4
        %p.r = [-R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.r = [-R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.rb = [-R*sin(p.U).*sin(p.V), R*sin(p.U).*cos(p.V), R*cos(p.U)];
        p.r_ = [-R*sin(p.u_).*sin(p.v_), R*sin(p.u_).*cos(p.v_), R*cos(p.u_)];
    end

    
p = p.update();    %updating x, y, z from r and pou
p = p.updateStale();
%p.J = abs(R^2*sin(p.u)); %Setting jacobian determinant

end


