function node = patchParameterise(r, p_in, p_out)
%Function to find the u-v coordinates of a point with cartesian coordinates r in p_in p_in
% in p_in p_out
% if r is not in p_in p_out returns (u, v) = (-1, -1).
    R = p_out.R;
    
    if p_out.numPatch==1
        r0 = [0,0,R];
        if r(2) < 0
            u = -1;
            v = -1;
        else 
            u = angle(r, r0);
            v = angle(r - norm(r)*cos(u)*[0,0,1], [1,0,0]);
        end
    elseif p_out.numPatch ==2
        r0 = [0,0,R];
        if r(2) > 0
            u = -1;
            v = -1;
        else
            u = angle(r0, r);
            v = angle(r - norm(r)*cos(u)*[0,0,1], [1,0,0]);
        end
    elseif p_out.numPatch==3
        r0 = [0,R,0];
        if r(3) < 0
            u = -1;
            v = -1;
        else
            u = angle(r0, r);
            v = angle(r - norm(r)*cos(u)*[0,1,0], [1,0,0]);
        end
    elseif p_out.numPatch==4
        r0 = [0,R,0];
        if r(3) > 0
            u = -1;
            v = -1;
        else
            u = angle(r0, r);
            v = angle(r-norm(r)*cos(u)*[0,1,0], [1,0,0]);
        end
    elseif p_out.numPatch==5
        r0 = [0,0,R];
        if(r(1) < 0)
            u = -1;
            v = -1;
        else
            u = angle(r0, r);
            v = angle(r - norm(r)*cos(u)*[0,0,1], [0,1,0]);
        end
    elseif p_out.numPatch==6
        r0 = [0,0,R];
        if(r(1) > 0)
           u = -1;
           v = -1;
        else
            u = angle(r0, r);
            v = angle(r - norm(r)*cos(u)*[0,0,1], [0,1,0]);
        end
        
        
    end
    
    if(u < 0 || u > pi)
        u = -1;
    end
 
    if(v < 0 || v > pi)
        v = -1;
    end
    node = [u, v];

end