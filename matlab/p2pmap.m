function node = p2pmap(r, p_in, p_out)
%Function to find the u-v coordinates of a point with cartesian coordinates r in p_in 
% in p_in p_out
% this includes points not in p_out, i.e v > pi is allowed 
    R = p_out.R;
    
    if p_out.numPatch==1
        r0 = [0,0,R];
        
        
            u = angle(r, r0);
            a = r - norm(r)*cos(u)*[0,0,1];
            b = [1, 0, 0];
            v = angle(a, b);
            if r(2) < 0
                %u = -inf;
                %v = -inf;
                v = 2*pi-v;
                
            end
%             if dot(cross(b, a), [0, 0, 1]) >= 0
%                 v = v;
%             else
%                 v = -v;
%             end
%                 
        
    elseif p_out.numPatch ==2
        r0 = [0,0,R];
        
            u = angle(r0, r);
            a = r - norm(r)*cos(u)*[0,0,1];
            b = [-1,0,0];
            v = angle(a, b);
            if r(2) > 0
                %u = -inf;
                %v = -inf;
                v = 2*pi-v;
            end
%             if dot(cross(b, a), [0, 0, 1]) >= 0
%                 v = v;
%             else
%                 v = -v;
%             end
        
    elseif p_out.numPatch==5
        r0 = [0,-R,0];
        
            u = angle(r0, r);
            a = r - norm(r)*cos(u)*[0,-1,0];
            b = [1,0,0];
            v = angle(a, b);
            if r(3) < 0
                %u = -inf;
                %v = -inf;
                v = 2*pi-v;
            end
%             if dot(cross(b, a), [0, 0, 1]) >= 0
%                 v = v;
%             else
%                 v = -v;
%             end        

    elseif p_out.numPatch==6
        r0 = [0,R,0];

        u = angle(r0, r);
        v = angle(r-norm(r)*cos(u)*[0,1,0], [1,0,0]);
        
        if r(3) > 0
            %u = -inf;
            %v = -inf;
            v = 2*pi-v;
        end
        
    elseif p_out.numPatch==3
        r0 = [0,0,R];

            u = angle(r0, r);
            a = r - norm(r)*cos(u)*[0,0,1];
            b = [0,-1,0];
            v = angle(a, b);
            if r(1) < 0
                %u = -inf;
                %v = -inf;
                v = 2*pi-v;
            end
%             if dot(cross(b, a), [0, 0, 1]) >= 0
%                 v = v;
%             else
%                 v = -v;
%             end            
%         
    elseif p_out.numPatch==4
        r0 = [0,0,R];
            u = angle(r0, r);
            a = r - norm(r)*cos(u)*[0,0,1];
            b = [0,1,0];
            v = angle(a, b);
            if r(1) > 0
                %u = -inf;
                %v = -inf;
                v = 2*pi-v;
            end
%             if dot(cross(b, a), [0, 0, 1]) >= 0
%                 v = v;
%             else
%                 v = -v;
%             end        
%         
        
    end
    
    if(v>pi)
        %u = -inf;
        %v = -inf;
    
%     elseif (v >= 0 & v <= pi)
%         u = u;
%         v = v;
%         
%     elseif (v >= -p_out.eps_strip & v < 0)
%         u = u;
%         v = v;
%         
%     elseif (v > -pi & v < -pi + p_out.eps_strip)
%         u = u;
%         v = 2*pi + v;
%         
%     else
%         u = -inf;
%         v = -inf
    end
    
    node = [u, v];

end