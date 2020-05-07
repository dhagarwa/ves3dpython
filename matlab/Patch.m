classdef Patch < handle
   properties
        u
        v
        nodes 
        numPatch
        numPatches
        neighbors
        numNodes
        Nu %# of u nodes
        Nv %# of v nodes
        h_u
        h_v
        x
        y
        z
        r
        J
        pou
        R
        q_dl %density double layer
        q_sl %density single layer
        n %unit normal
        eps_strip %extend hemispherical patch by eps_strip more v coordinate
        linkPatch
        links %u-v parameters of points in current patch in other patch
   end
   
   methods
      function obj = Patch(m, n, R, numPatch, neighbors)
         obj.Nu = m;
         obj.Nv = n;
         obj.eps_strip = pi/6;
         [obj.u,obj.v]=ndgrid(pi*(1:m)/(m+1), -obj.eps_strip + (pi + 2*obj.eps_strip)*(1:n)/(n+1)); %u and v both are m x n matrices
         obj.u = obj.u(:);
         obj.v = obj.v(:);
         obj.nodes =[obj.u,obj.v]; % N by 2 matrix listing x,y coordinates of all N=m*n nodespi
         obj.h_u = pi/ (m+1);
         obj.h_v = (pi + 2*obj.eps_strip)/ (n+1);    
         obj.numNodes = m*n; %no of nodes on one patch in u-v coordinates
         obj.numPatch = numPatch;
         obj.neighbors = neighbors;
         obj.R = R;
         obj.numPatches = 4; 
         obj.links = -1*ones(obj.numNodes, 2*obj.numPatches);
      end
      
%       function initialiseMaps(obj, param, jac)
%           obj.jacobian = jac;
%           obj.param = param;
%           
%       end
      
      function TotalNodes = getTotalNodes(objArray)
         TotalNodes = 0;
         for i = 1:numel(objArray)
            TotalNodes = TotalNodes + objArray(i).numNodes;
         end
      end
      
      function obj = update(obj)
         obj.x = obj.r(:, 1);
         obj.y = obj.r(:, 2);
         obj.z = obj.r(:, 3);
         obj.n = obj.r/obj.R;
         qk = zeros(obj.numNodes, 1);
         for pp= 1:obj.numPatches
             qk = qk+ obj.spou(pp, obj.r);
         end
         obj.pou = obj.spou(obj.numPatch, obj.r)./qk;
          
      end

    function val = spou(obj, patch, r) %patch is patch num here
    [v1, u1, rho1] = cart2sph(r(:,1), r(:,2), r(:,3));
    %u = pi/2 - u;
    %u; v;
    if patch==1
        r0 = [0,obj.R,0];
    elseif patch ==2    
        r0 = [0,-obj.R,0];
    elseif patch==5
        r0 = [0,0,obj.R];
    elseif patch==6
        r0 = [0,0,-obj.R];
    elseif patch==3
        r0 = [obj.R,0,0];
    elseif patch==4
        r0 = [-obj.R,0,0];
        
    end
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    %u0 = pi/2-u0;
    %u0; v0;
    d = (6/5) *(5/12) * pi*obj.R;
    nn = size(r, 1);
    val = zeros(nn, 1);
    for ii=1:nn
        t = greatCircleDistance(u0, v0, u1(ii), v1(ii), obj.R)/d;
        if t>=1
            val(ii) = 0;
        elseif t ==0
            val(ii) = 1;
        else
            val(ii) = exp((2*exp(-1/t))/ (t-1));
        end
    end
    
    
    
    end
    
    function [f_du, f_dv] = grad(obj, f)
         
%         f_mat = reshape(f, [obj.Nu obj.Nv]); %u varies column wise, v rowwise
%         f_mat = [zeros(1, obj.Nv);f_mat];
%         f_mat = [zeros(obj.Nu+1, 1) f_mat];
%         f_mat = 2*f_mat; %2 for scaling due to transformation [0, pi]->[-pi, pi]
        f = f.*obj.pou; %partition of unity done
        f_mat = ves2fft(f, obj.Nu, obj.Nv);
        [f_du, f_dv] = deriv(f_mat); %transpose because fft code has u vlues in row and v in column
         f_du = fft2ves(f_du, obj.Nu, obj.Nv);
         f_dv = fft2ves(f_dv, obj.Nu, obj.Nv);
        
    end

        function [f_du, f_dv] = grad_cheb(obj, f)
         
%         f_mat = reshape(f, [obj.Nu obj.Nv]); %u varies column wise, v rowwise
%         f_mat = [zeros(1, obj.Nv);f_mat];
%         f_mat = [zeros(obj.Nu+1, 1) f_mat];
%         f_mat = 2*f_mat; %2 for scaling due to transformation [0, pi]->[-pi, pi]
        
        f_mat = ves2fft(f, obj.Nu, obj.Nv);
        [f_du, f_dv] = deriv(f_mat); %transpose because fft code has u vlues in row and v in column
         f_du = fft2ves(f_du, obj.Nu, obj.Nv);
         f_dv = fft2ves(f_dv, obj.Nu, obj.Nv);
        
        end
    
       function [f_du, f_dv] = interpolate(obj, f)
         
        f = f.*obj.pou; %partition of unity done
        f_mat = ves2fft(f, obj.Nu, obj.Nv);
        [f_du, f_dv] = deriv(f_mat); %transpose because fft code has u vlues in row and v in column
         f_du = fft2ves(f_du, obj.Nu, obj.Nv);
         f_dv = fft2ves(f_dv, obj.Nu, obj.Nv);
        
    end
    
   end
end


