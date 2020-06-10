classdef Patch < handle
   properties
        u
        v
        nodes
        U
        V
        Nodes
        numPatch
        numPatches
        neighbors
        numNodes
        Nu %# of u nodes
        Nv %# of v nodes
        NU % # of v nodes
        NV % # of U nodes
        h_u
        h_v
        x
        y
        z
        X %x with boundary
        Y
        Z
        Xr %reference x
        Yr % " y
        Zr % " z
        r
        rb % r but  with boundary nodes included
        rrb % r in reference configuration with boundary
        J  % current jacobian
        Jr %reference configuration jacobian
        pou
        R
        q_dl %density double layer
        q_sl %density single layer
        n %unit normal current
        nr %unit normal reference configuration
        eps_strip %extend hemispherical patch by eps_strip more v coordinate
        linkPatch
        links %u-v parameters of points in current patch in other patch
        Du % u-derivative FDM matrix 
        Dv % v-derivative FDM matrix
        ooa %order of accuracy for FDM
        %Reference tangents
        Xr_u %boundary included
        Xr_v
        Yr_u
        Yr_v
        Zr_u
        Zr_v  
        rr_u
        rr_v
        %Tangents
        X_u %boundary included
        X_v
        Y_u
        Y_v
        Z_u
        Z_v
        r_u
        r_v
        %First fundamental form
        E
        F
        G
        W
        
   end
   
   methods
      function obj = Patch(m, n, R, numPatch, neighbors)
         obj.Nu = m;
         obj.Nv = n;
         obj.eps_strip = 0;
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
         obj.numPatches = 6; 
         obj.links = -1*ones(obj.numNodes, 2*obj.numPatches);
         %FDM full U, V including boundary of u-v patch
         [obj.U,obj.V]=ndgrid(pi*(1:m)/(m+1), -obj.eps_strip + (pi + 2*obj.eps_strip)*(1:n)/(n+1)); %u and v both are m x n matrices
         obj.U = obj.U(:);
         obj.V = obj.V(:);
         obj.Nodes =[obj.U,obj.V]; % N by 2 matrix listing x,y coordinates of all N=m*n nodespi
         obj.NU = m;
         obj.NV = n;
         obj.ooa = 8;
         [obj.Du, obj.Dv] = getNonCompactFDmatrix2D(obj.NU, obj.NV, obj.h_u, obj.h_v, 1, obj.ooa);
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
         obj.X = obj.rb(:, 1);
         obj.Y = obj.rb(:, 2);
         obj.Z = obj.rb(:, 3); 
         obj.rrb = obj.rb;
         obj.Xr = obj.rrb(:, 1);
         obj.Yr = obj.rrb(:, 2);
         obj.Zr = obj.rrb(:, 3);
         [obj.Xr_u, obj.Xr_v ] = obj.grad_FDM(obj.Xr);
         [obj.Yr_u, obj.Yr_v ] = obj.grad_FDM(obj.Yr);
         [obj.Zr_u, obj.Zr_v ] = obj.grad_FDM(obj.Zr);
         obj.rr_u = [obj.Xr_u, obj.Yr_u, obj.Zr_u];

         obj.rr_v = [obj.Xr_v, obj.Yr_v, obj.Zr_v];
         
         %set new jacobian determinant J ||r_u \times r_v||
         
         %set normal
         temp_cross = cross(obj.rr_u, obj.rr_v, 2);
         obj.Jr = abs(vecnorm(temp_cross, 2, 2));
         obj.nr = diag(vecnorm(temp_cross, 2, 2).^(-1))*temp_cross;
         qk = zeros(obj.numNodes, 1);
         for pp= 1:obj.numPatches
             qk = qk+ obj.spou(pp, obj.r);
         end
         obj.pou = obj.spou(obj.numPatch, obj.r)./qk; %TODO: obj.r is incorrect if initial shape is non sphere
%          plot(obj.spou(obj.numPatch, obj.r));
%          hold on
%          plot(qk)
%          hold on
%          plot(obj.pou)
%          legend("spou", 'qk', 'pou');
%          disp("pou component");
      end

       function obj = updateStale(obj) %update to be called after patch.rb changes
         obj.x = obj.r(:, 1); %new positions
         obj.y = obj.r(:, 2);
         obj.z = obj.r(:, 3);
         obj.X = obj.rb(:, 1); %new positions with boundary
         obj.Y = obj.rb(:, 2);
         obj.Z = obj.rb(:, 3); 
         %set new jacobian i.e. tangents x_u, y_u, z_u
         [obj.X_u, obj.X_v ] = obj.grad_FDM(obj.X);
         [obj.Y_u, obj.Y_v ] = obj.grad_FDM(obj.Y);
         [obj.Z_u, obj.Z_v ] = obj.grad_FDM(obj.Z);
         obj.r_u = [obj.X_u, obj.Y_u, obj.Z_u];

         obj.r_v = [obj.X_v, obj.Y_v, obj.Z_v];
         
         %set new jacobian determinant J ||r_u \times r_v||
         
         %set normal
         temp_cross = cross(obj.r_u, obj.r_v, 2);
         obj.J = abs(vecnorm(temp_cross, 2, 2));
         obj.n = diag(vecnorm(temp_cross, 2, 2).^(-1))*temp_cross;
         [obj.E, obj.F, obj.G, obj.W] = obj.ffform();
          
         
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
    d = (5/12) * pi*obj.R;
    nn = size(r, 1);
    val = zeros(nn, 1);
    t_vec = zeros(nn, 1);
    for ii=1:nn
        t = greatCircleDistance(u0, v0, u1(ii), v1(ii), obj.R)/d;
        t_vec(ii) = t;
        if t>=1
            val(ii) = 0;
        elseif t ==0
            val(ii) = 1;
        else
            val(ii) = exp((2*exp(-1/t))/ (t-1));
        end
    end
    
    %plot(val)
    %disp("plotting great circle distance");
    end
    
    function [f_du, f_dv] = grad(obj, f)
        [fphi_du, fphi_dv] = obj.grad_patch(f);
        fphi_du;
        [phi_du, phi_dv] = obj.grad_patch(ones(size(obj.numNodes, 1), 1));
        phi_du;
        phi = obj.pou;
        phi(phi ==0) = 1; %To avoid division by zero
        f_du = (fphi_du - f.*phi_du)./obj.pou ; 
        f_dv = (fphi_dv - f.*phi_dv)./obj.pou ;
        
    end
    
    function [fphi_du, fphi_dv] = grad_patch(obj, f)
         
%         f_mat = reshape(f, [obj.Nu obj.Nv]); %u varies column wise, v rowwise
%         f_mat = [zeros(1, obj.Nv);f_mat];
%         f_mat = [zeros(obj.Nu+1, 1) f_mat];
%         f_mat = 2*f_mat; %2 for scaling due to transformation [0, pi]->[-pi, pi]
        f = f.*obj.pou; %partition of unity done
        %plot(f)
        f_mat = ves2fft(f, obj.Nu, obj.Nv);
        %f_mat(15, 2:(obj.Nu+1))
        
         %plot(f_mat(15, 2:(obj.Nv+1)))
         %hold on
        [fphi_du, fphi_dv] = deriv(f_mat); % fft code has u values in row and v in column
        
         fphi_du = 2 * fft2ves(fphi_du, obj.Nu +1, obj.Nv +1); %2 is scaling factor =  2*pi/ pi
         %plot(f_du((14*obj.Nu + 1):(15*obj.Nu)))
         
         %plot(f_du - 2*cos(2*obj.u));
         %hold on
         fphi_dv = (2*pi/(pi + 2*obj.eps_strip)) * fft2ves(fphi_dv, obj.Nu+1, obj.Nv+1); %(2*pi/(pi + obj.eps_strip)) is scaling factor
         %f_dv
         %plot(f_dv  - 2*pi/(2*obj.eps_strip + pi)*cos(2*pi/(2*obj.eps_strip + pi)*(obj.v + obj.eps_strip)));
         %f_dv;
         
    end


    
    %Finite difference derivative
     function [f_du, f_dv] = grad_FDM(obj, f) %ood = order of derivative, ooa = order of accuracy
         
        f_fdm = ves2fdm(f, obj.NU, obj.NV); %transpose because fdm code has U values in row and V in column
        %[f_du, f_dv] = derivFDM(f_fdm, obj.NU, obj.NV, obj.h_u, obj.h_v, ood, ooa); 
        f_du = obj.Du*f_fdm;
        f_dv = obj.Dv*f_fdm;
        f_du = fdm2ves(f_du, obj.NU, obj.NV);
        f_dv = fdm2ves(f_dv, obj.NU, obj.NV);
        
     end
        
       %Incomplete function
    function [f_du, f_dv] = interpolate(obj, f)
         
        f = f.*obj.pou; %partition of unity done
        f_mat = ves2fft(f, obj.Nu, obj.Nv);
        [f_du, f_dv] = deriv(f_mat); %transpose because fft code has u values in row and v in column
         f_du = fft2ves(f_du, obj.Nu+1, obj.Nv+1);
         f_dv = fft2ves(f_dv, obj.Nu+1, obj.Nv+1);
        
    end
    
    %Function to calculate first fundamental form
    function [E, F, G, W] = ffform(obj)
       E = dot(obj.r_u, obj.r_u, 2); 
       F = dot(obj.r_u, obj.r_v, 2); 
       G = dot(obj.r_v, obj.r_v, 2);
       W = (E.*G - F.^2).^(0.5);
    end
    
    function val = surf_div(obj, f)
      %Function to calculate surface divergence
      %of vector valued  f. size(f) = numNodes x 3
        [fx_du, fx_dv] = obj.grad_FDM(f(:, 1));
        [fy_du, fy_dv] = obj.grad_FDM(f(:, 2));
        [fz_du, fz_dv] = obj.grad_FDM(f(:, 3));
        
        f_du = [fx_du, fy_du, fz_du];
        f_dv = [fx_dv, fy_dv, fz_dv];
        Gf_du = [obj.G./(obj.W).^2 obj.G./(obj.W).^2 obj.G./(obj.W).^2].*f_du;
        Ff_dv = [obj.F./(obj.W).^2 obj.F./(obj.W).^2 obj.F./(obj.W).^2].*f_dv;
        Ef_dv = [obj.E./(obj.W).^2 obj.E./(obj.W).^2 obj.E./(obj.W).^2].*f_dv;
        Ff_du = [obj.F./(obj.W).^2 obj.F./(obj.W).^2 obj.F./(obj.W).^2].*f_du;
        a = Gf_du - Ff_dv;
        b = Ef_dv - Ff_du;
        val = dot(a, obj.r_u, 2) + dot(b, obj.r_v, 2);
    end
    
    %Function to calculate shear force
    function fs = shearForce(obj, Es, Ed) %Es: shear modulus, Ed: dilatation modulus
        %given by divergence of symmetric part of stress tensor
%          r_u = [obj.X_u, obj.Y_u, obj.Z_u];
%          r_v = [obj.X_v, obj.Y_v, obj.Z_v];
%          rr_u = [obj.Xr_u, obj.Yr_u, obj.Zr_u];
%          rr_v = [obj.Xr_v, obj.Yr_v, obj.Zr_v];
         
         F = []; %relative surface deformation gradient
         V = []; %left Cauchy Green deformation tensor = F*F'
         
         P = []; % I - n*n', the projection on tangent plane
         first_invar = []; %I1 
         sec_invar = []; % I2
         tau = [];
       for ii=1:obj.numNodes
           
          a_r = [obj.rr_u(ii,:)' obj.rr_v(ii,:)' obj.nr(ii,:)']; %reference tangents
          a_c = [obj.r_u(ii,:)' obj.r_v(ii,:)' 0*obj.n(ii,:)']; %current of infinitesimal elements
          F_ii = (a_r'\a_c')';
          V_ii = F_ii*F_ii';
          [EV_ii, D_ii] = eig(V_ii);
          lambda_sq_ii = sort(diag(D_ii));
          lambda1_sq_ii = lambda_sq_ii(3);
          lambda2_sq_ii = lambda_sq_ii(2);
          first_invar_ii = lambda1_sq_ii + lambda2_sq_ii - 2;
          sec_invar_ii = lambda1_sq_ii * lambda2_sq_ii - 1;
          P_ii = eye(3) - obj.n(ii, :)'*obj.n(ii, :);
          Js_ii = sqrt(lambda1_sq_ii * lambda2_sq_ii);
          tau_ii = Es/(2*Js_ii)*(first_invar_ii + 1) * V_ii + (Js_ii/2)*(Ed*sec_invar_ii - Es)*P_ii; %symmetric stress tensor
          %norm(F_ii*a_r - a_c)
          F = [F; F_ii];
          V = [V; V_ii];
          P = [P; P_ii];
          first_invar = [first_invar; first_invar_ii];
          sec_invar = [sec_invar; sec_invar_ii];
          tau = [tau; tau_ii];
       end
       nrows = size(tau, 1);
       taux = tau(1:3:nrows, :);
       tauy = tau(2:3:nrows, :);
       tauz = tau(3:3:nrows, :);
       fsx = obj.surf_div(taux);
       fsy = obj.surf_div(tauy);
       fsz = obj.surf_div(tauz);
       fs = [fsx, fsy, fsz]; 
    end
   end
end


