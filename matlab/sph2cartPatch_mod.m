function val = sph2cartPatch_mod(u, v, numPatch)
    
    R = 1;

    %u(u == -1) = NaN;
    %v(v == -1) = NaN;
    if numPatch==1 
        val = [R*sin(u).*cos(v), R*sin(u).*sin(v), R*cos(u)];
    elseif numPatch==2
        val = [-R*sin(u).*cos(v), -R*sin(u).*sin(v), R*cos(u)];
    elseif numPatch==5
        val = [R*sin(u).*cos(v), -R*cos(u), R*sin(u).*sin(v)];

    elseif numPatch==6
        val = [R*sin(u).*cos(v), R*cos(u), -R*sin(u).*sin(v)];

    elseif numPatch==3
        val = [R*sin(u).*sin(v), -R*sin(u).*cos(v), R*cos(u)];

    elseif numPatch==4
        val = [-R*sin(u).*sin(v), R*sin(u).*cos(v), R*cos(u)];
    end



end