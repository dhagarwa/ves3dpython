function p = standardSpherePatch(m, n, numPatch, R)


p = Patch(m, n, R, numPatch, []);

    if numPatch==1 
        p.r = [R*sin(p.u).*cos(p.v), R*sin(p.u).*sin(p.v), R*cos(p.u)];
    elseif numPatch==2
        p.r = [R*sin(p.u).*cos(p.v), -R*sin(p.u).*sin(p.v), R*cos(p.u)];
    elseif numPatch==3
        p.r = [R*sin(p.u).*cos(p.v), R*cos(p.u), R*sin(p.u).*sin(p.v)];

    elseif numPatch==4
        p.r = [R*sin(p.u).*cos(p.v), R*cos(p.u), -R*sin(p.u).*sin(p.v)];

    elseif numPatch==5
        p.r = [R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];


    elseif numPatch==6
        p.r = [-R*sin(p.u).*sin(p.v), R*sin(p.u).*cos(p.v), R*cos(p.u)];
    end

    
p.update();    %updating x, y, z from r
p.J = abs(R^2*sin(p.u)); %Setting jacobian


end


