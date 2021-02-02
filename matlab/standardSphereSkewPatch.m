function p = standardSphereSkewPatch(m, n, numPatch, R, skew)


p = Patch(m, n, R, numPatch, []);

    if numPatch==1 
        
        u = pi*1./(1 + (p.u./(pi-p.u)).^(-skew));
        u_ = pi*1./(1 + (p.u_./(pi-p.u_)).^(-skew));
        U = pi*1./(1 + (p.U./(pi-p.U)).^(-skew));
        p.r = [R*sin(u).*cos(p.v), R*sin(u).*sin(p.v), R*cos(u)];
        p.r_ = [R*sin(u_).*cos(p.v_), R*sin(u_).*sin(p.v_), R*cos(u_)];
        p.rb = [R*sin(U).*cos(p.V), R*sin(U).*sin(p.V), R*cos(U)];
    elseif numPatch==2
        %p.r = [R*sin(p.u).*cos(p.v), -R*sin(p.u).*sin(p.v), R*cos(p.u)];
        %before patch u-v direction alignment 
        u = pi*1./(1 + (p.u./(pi-p.u)).^(-skew)); 
        u_ = pi*1./(1 + (p.u_./(pi-p.u_)).^(-skew));
        U = pi*1./(1 + (p.U./(pi-p.U)).^(-skew));        
        p.r = [-R*sin(u).*cos(p.v), -R*sin(u).*sin(p.v), R*cos(u)];
        p.r_ = [-R*sin(u_).*cos(p.v_), -R*sin(u_).*sin(p.v_), R*cos(u_)];
        p.rb = [-R*sin(U).*cos(p.V), -R*sin(U).*sin(p.V), R*cos(U)];
    elseif numPatch==5
        u = pi*1./(1 + (p.u./(pi-p.u)).^(-skew)); 
        u_ = pi*1./(1 + (p.u_./(pi-p.u_)).^(-skew));
        U = pi*1./(1 + (p.U./(pi-p.U)).^(-skew));        
        p.r = [R*sin(u).*cos(p.v), -R*cos(u), R*sin(u).*sin(p.v)]; %u 0 to pi, v 0 to pi
        %p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)]; % u 0 to pi/2 , v 0 to 2*pi, patch is to be extended by mirror image
        p.rb = [R*sin(U).*cos(p.V), -R*cos(U), R*sin(U).*sin(p.V)];
        p.r_ = [R*sin(u_).*cos(p.v_), -R*cos(u_), R*sin(u_).*sin(p.v_)];
    elseif numPatch==6
        u = pi*1./(1 + (p.u./(pi-p.u)).^(-skew)); 
        u_ = pi*1./(1 + (p.u_./(pi-p.u_)).^(-skew));
        U = pi*1./(1 + (p.U./(pi-p.U)).^(-skew));         
        p.r = [R*sin(u).*cos(p.v), R*cos(u), -R*sin(u).*sin(p.v)]; %u 0 to pi, v 0 to pi
        %p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)]; % u pi/2 to pi, v 0 to 2*pi, patch is to be extended by mirror image
        p.rb = [R*sin(U).*cos(p.V), R*cos(U), -R*sin(U).*sin(p.V)];
        p.r_ = [R*sin(u_).*cos(p.v_), R*cos(u_), -R*sin(u_).*sin(p.v_)];
    elseif numPatch==3
        u = pi*1./(1 + (p.u./(pi-p.u)).^(-skew)); 
        u_ = pi*1./(1 + (p.u_./(pi-p.u_)).^(-skew));
        U = pi*1./(1 + (p.U./(pi-p.U)).^(-skew));         
        %p.r = [R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.r = [R*sin(u).*sin(p.v), -R*sin(u).*cos(p.v), R*cos(u)];
        p.rb = [R*sin(U).*sin(p.V), -R*sin(U).*cos(p.V), R*cos(U)];
        p.r_ = [R*sin(u_).*sin(p.v_), -R*sin(u_).*cos(p.v_), R*cos(u_)];
    elseif numPatch==4
        u = pi*1./(1 + (p.u./(pi-p.u)).^(-skew)); 
        u_ = pi*1./(1 + (p.u_./(pi-p.u_)).^(-skew));
        U = pi*1./(1 + (p.U./(pi-p.U)).^(-skew));         
        %p.r = [-R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
        p.r = [-R*sin(u).*sin(p.v), R*sin(u).*cos(p.v), R*cos(u)];
        p.rb = [-R*sin(U).*sin(p.V), R*sin(U).*cos(p.V), R*cos(U)];
        p.r_ = [-R*sin(u_).*sin(p.v_), R*sin(u_).*cos(p.v_), R*cos(u_)];
    end

    
p = p.update();    %updating x, y, z from r and pou
p = p.updateStale();
%p.J = abs(R^2*sin(p.u)); %Setting jacobian determinant

end


