classdef Surface < handle
   properties
        patches
        link %1 linked , 0 not linked
        link_
        p2plinked
        derivlinked 
        secondderivlinked
        poulinked
        numPatches
        kb %bending modulus
        Es %shear modulus
        Ed %dilatation modulus
   end
   
   methods
      function obj = Surface(patches, vesprops)
            obj.patches = patches;
            obj.numPatches = length(patches);
            obj.link = obj.linkPatches();
            obj.link_ = obj.linkPatches_();
            obj.p2plinked = obj.linkp2p(); 
            obj.derivlinked = obj.derivp2p(); %need to chamge for Symm patch
            %obj.secondderivlinked = obj.secondderivp2p(); %need to chamge for Symm patch
            obj.kb = vesprops(1);
            obj.Es = vesprops(2);
            obj.Ed = vesprops(3);
            
      end
      
      function vals = getPosition(obj)
          vals = [];
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              vals = [vals;patch.r];              
          end
          
      end
      
      function [] = interfacialForce(obj)
          %setting interfacial force in each patch
          m = obj.patches(1).Nu; n = obj.patches(1).Nv;
          f = [];
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              f_ii = patch.shearForce(obj.Es, obj.Ed);
              f = [f; f_ii];
              %patch.q_sl = f_ii;
              %disp('Force');
          end
          
          fx_blend = obj.blendSurfaceFunction(f(:,1));
          fy_blend = obj.blendSurfaceFunction(f(:,2));
          fz_blend = obj.blendSurfaceFunction(f(:,3));
          
          f_blend = [fx_blend fy_blend fz_blend];
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              f_ii = f_blend((ii-1)*m*n+1:ii*m*n, :);
              patch.q_sl = f_ii;
              %disp('Force');
          end          
      end
      
      function [] = updateSurface(obj, u, dt)          
          %updating stale patches
          blk_size = size(u, 1)/obj.numPatches;
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              patch.r = patch.r + u((ii-1)*blk_size + 1:ii*blk_size, :)*dt;
              patch.rb = patch.r;
                             
          end     
          
      end
            
      
      function [] = updateStale(obj)          
          %updating stale patches
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              patch = patch.updateStale();
               
          end     
          
      end
      
      function center = getCenter(obj)
          r = obj.getPosition();
          center = sum(r)/size(r, 1);
          
      end

      function x_disp = recenter_x(obj) %recenter x to 0  
          center = obj.getCenter();
          x_disp = center(1)-2;
          disp = [ones(size(obj.patches(1).r, 1), 1) zeros(size(obj.patches(1).r, 1), 1) zeros(size(obj.patches(1).r, 1), 1)];
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);  
              %patch.r
              %x_disp*disp
              patch.r = patch.r - x_disp*disp;            
              patch.rb = patch.r;
                             
          end              
          
      end      
      
      
      function [] = blendSurface(obj)
          
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
              x = zeros(patch.numNodes,1);
              y = zeros(patch.numNodes,1);
              z = zeros(patch.numNodes,1);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    x_linked_patch = linked_patch.x.*linked_patch.pou;
                    y_linked_patch = linked_patch.y.*linked_patch.pou;
                    z_linked_patch = linked_patch.z.*linked_patch.pou;
                    
                    xq = patch.links(:, 2*jj-1:2*jj); %query points from links
                    x_linked_interp = interpPatch(x_linked_patch, xq, linked_patch);
                    y_linked_interp = interpPatch(y_linked_patch, xq, linked_patch);
                    z_linked_interp = interpPatch(z_linked_patch, xq, linked_patch);
                    
                    x = x + x_linked_interp;
                    y = y + y_linked_interp;
                    z = z + z_linked_interp;
                     
                    
               end
            end
            x = patch.x.*patch.pou + x;
            y = patch.y.*patch.pou + y;
            z = patch.z.*patch.pou + z;
            patch.r = [x y z];
            patch.rb = patch.r;

          end          
      end
      
      function val = linkPatches(obj)
         for ii=1:length(obj.patches)
            patch = obj.patches(ii);
            for jj=1:patch.numNodes
               r0 = patch.sphcart(jj, :);
               for kk=1:length(obj.patches)
                  p = obj.patches(kk);
                  if kk ~= patch.numPatch
                     %find out u-v coordinates in patch p 
                     
                      node = patchParameteriseSymm(r0, patch, p);
                      patch.links(jj, 2*(kk-1) + 1:2*kk) = node;
                  end
               end
            end
         end
         %obj.patches(1).links
          val = 1;
      end

      function val = linkPatches_(obj)
         for ii=1:length(obj.patches)
            patch = obj.patches(ii);
            for jj=1:patch.numNodes_
               r0 = patch.sphcart_(jj, :);
               for kk=1:length(obj.patches)
                  p = obj.patches(kk);
                  if kk ~= patch.numPatch
                     %find out u-v coordinates in patch p 
                     
                      node = patchParameteriseSymm(r0, patch, p);
                      patch.links_(jj, 2*(kk-1) + 1:2*kk) = node;
                  end
               end
            end
         end
         %obj.patches(1).links
          val = 1;
      end
      
      
      %function to get u-v coordinates of nodes in one patch according to
      %other patch without -inf if not present in patch
      function val = linkp2p(obj)
         for ii=1:length(obj.patches)
            patch = obj.patches(ii);
            for jj=1:patch.numNodes
               r0 = patch.sphcart(jj, :);
               for kk=1:length(obj.patches)
                  p = obj.patches(kk);
                  
                  %find out u-v coordinates in patch p 
                  node = p2pmapSymm(r0, patch, p);
                  patch.cl_links(jj, 2*(kk-1) + 1:2*kk) = node;
                  
               end
            end
         end
         %obj.patches(1).links
          val = 1;
      end

     
 
      %function to get u-v coordinates of nodes in one patch according to
      %other patch without -inf if not present in patch
      function val = derivp2p(obj)
         for ii=1:length(obj.patches)
            p1 = obj.patches(ii);
            for jj=1:length(obj.patches)
                p2 = obj.patches(jj);
                  
                %find out u-v deriv  in patch p2 wrt p1 u-v 
                pc_deriv = p2pderivmap(p1, p2);
                p1.deriv_links(:, 4*(jj-1) + 1:4*jj) = pc_deriv;
                pc_deriv = p2pderivmap_(p1, p2);
                p1.deriv_links_(:, 4*(jj-1) + 1:4*jj) = pc_deriv;
                         
               
            end
         end
         %obj.patches(1).links
          val = 1;
      end  
      
      function val = secondderivp2p(obj)
         for ii=1:length(obj.patches)
            p1 = obj.patches(ii);
            for jj=1:length(obj.patches)
                p2 = obj.patches(jj);
                  
                %find out second u-v deriv  in patch p2 wrt p1 u-v 
                pc_secondderiv = p2psecondderivmap(p1, p2);
                p1.secondderiv_links(:, 8*(jj-1) + 1:8*jj) = pc_secondderiv;
                %pc_secondderiv = p2psecondderivmap_(p1, p2);
                %p1.secondderiv_links_(:, 4*(jj-1) + 1:4*jj) = pc_deriv;
                         
               
            end
         end
         %obj.patches(1).links
          val = 1;
      end  
      
      
      function [f_du, f_dv] = deriv(obj, f)
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          f_du = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dv = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_du_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dv_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
          for ii=1:length(obj.patches)
              patch = obj.patches(ii);
              f_patch = f(:, patch.numPatch);
              [f_du_patch, f_dv_patch] = patch.grad(f_patch);
              % Add contribution from other patches
              %TODO
  
              f_du_patches(:, patch.numPatch) = f_du_patch;
              f_dv_patches(:, patch.numPatch) = f_dv_patch;
              f_du(:, patch.numPatch) = f_du_patch;
              f_dv(:, patch.numPatch) = f_dv_patch;
 
              
          end
          
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    f_du_linked_patch = f_du_patches(:, jj);
                    f_dv_linked_patch = f_dv_patches(:, jj);
                    xq = patch.links(:, 2*jj-1:2*jj); %query points from links
                    f_du_linked_interp = interpPatch(f_du_linked_patch, xq, linked_patch);
                    f_dv_linked_interp = interpPatch(f_dv_linked_patch, xq, linked_patch);
                    f_du(:, ii) = f_du(:, ii) + f_du_linked_interp;
                    f_dv(:, ii) = f_dv(:, ii) + f_dv_linked_interp;  
                    
               end
            end
            

          end
      end

      
      function [f_du, f_dv] = derivFFT(obj, f)
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          f_du = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dv = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_du_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dv_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
           pn = 5;
          for ii=1:length(obj.patches)
              patch = obj.patches(ii);
              f_patch = f(:, patch.numPatch);
              fphi_patch = f_patch.*patch.pou;
              if ii==pn
                 stem3(reshape(fphi_patch,[patch.Nu,patch.Nv])); 
              end
              [fphi_du_patch, fphi_dv_patch] = patch.grad_patch(fphi_patch);
              % Add contribution from other patches
              %TODO
  
              f_du_patches(:, patch.numPatch) = fphi_du_patch;
              f_dv_patches(:, patch.numPatch) = fphi_dv_patch;
              f_du(:, patch.numPatch) = fphi_du_patch;
              f_dv(:, patch.numPatch) = fphi_dv_patch;
 
              
          end
%           [fphi_du_fdm, fphi_dv_fdm] = patch.grad_FDM(fphi_patch);
%           error = max(abs(fphi_du_fdm - fphi_du_patch))
            
           figure;
           stem3(reshape(f(:,pn),[patch.Nu,patch.Nv]));
           figure;
           stem3(reshape(f_du(:,pn), [patch.Nu,patch.Nv]));
           figure;
           stem3(reshape(f_dv(:,pn), [patch.Nu,patch.Nv]));
%           figure;
%           stem3(reshape(fphi_dv_patch, [patch.Nu,patch.Nv]));
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    f_du_linked_patch = f_du_patches(:, jj);
                    f_dv_linked_patch = f_dv_patches(:, jj);
                    xq = patch.links(:, 2*jj-1:2*jj); %query points from links
                    u_linked_du = patch.deriv_links(:, 4*jj-3); %query points derivatives from links
                    u_linked_dv = patch.deriv_links(:, 4*jj-2); %query points derivatives from links
                    v_linked_du = patch.deriv_links(:, 4*jj-1); %query points derivatives from links
                    v_linked_dv = patch.deriv_links(:, 4*jj); %query points derivatives from links
                    %vq = patch.deriv_links(:, 2*jj); %query points derivatives from links
                    f_du_linked_interp = interpPatch2(f_du_linked_patch, xq, linked_patch);
                    f_dv_linked_interp = interpPatch2(f_dv_linked_patch, xq, linked_patch);
                    %[u_linked_du, u_linked_dv] = patch.grad_FDM(uq);
                    %[v_linked_du, v_linked_dv] = patch.grad_FDM(vq);
                    f_du(:, ii) = f_du(:, ii) + f_du_linked_interp.*u_linked_du + f_dv_linked_interp.*v_linked_du;
                    f_dv(:, ii) = f_dv(:, ii) + f_du_linked_interp.*u_linked_dv + f_dv_linked_interp.*v_linked_dv;  
                    
               end
            end
            

          end
      end
      
 
      function [f_du, f_dv] = derivFDM(obj, f)
          %Doesn't work for boundary of patch 
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          m = obj.patches(1).Nu; n = obj.patches(1).Nv;
          f_du = zeros((m+2)*(n+2), length(obj.patches));
          f_dv = zeros((m+2)*(n+2), length(obj.patches));
          f_du_patches = zeros((m+2)*(n+2), length(obj.patches));
          f_dv_patches = zeros((m+2)*(n+2), length(obj.patches));
           pn = 5;
          for ii=1:length(obj.patches)
              patch = obj.patches(ii);
              f_patch = f(:, patch.numPatch);
              %fphi_patch = f_patch.*patch.pou;
%               if ii==pn
%                  stem3(reshape(f_patch,[m+2,n+2])); 
%               end
              [f_du_patch, f_dv_patch] = patch.grad_FDM_(f_patch);
              % Add contribution from other patches
              %TODO
  
              f_du_patches(:, patch.numPatch) = f_du_patch;
              f_dv_patches(:, patch.numPatch) = f_dv_patch;
                
              f_du(:, patch.numPatch) = f_du_patch.*patch.pou_;
              f_dv(:, patch.numPatch) = f_dv_patch.*patch.pou_;
 
              if ii==pn                  
                 test_du = stripboundary(f_du_patch,patch); 
                 test_dv = stripboundary(f_dv_patch,patch);
                 test_du(12)
                 test_dv(12)
              end
              
              
          end
%           [fphi_du_fdm, fphi_dv_fdm] = patch.grad_FDM(fphi_patch);
%           error = max(abs(fphi_du_fdm - fphi_du_patch))
            
%             figure;
%             stem3(reshape(f_dv_patches(:,pn),[patch.Nu+2,patch.Nv+2]));
%             figure;
%             stem3(reshape(f_du(:,pn), [patch.Nu,patch.Nv]));
%             figure;
%             stem3(reshape(f_dv(:,pn), [patch.Nu,patch.Nv]));
%           figure;
%           stem3(reshape(fphi_dv_patch, [patch.Nu,patch.Nv]));
%             isnan(f_du_patches)
%             isnan(f_dv_patches)
%             isnan(f_du)
%             isnan(f_dv)
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    f_du_linked_patch = f_du_patches(:, jj);
                    f_dv_linked_patch = f_dv_patches(:, jj);
                    xq = patch.links_(:, 2*jj-1:2*jj); %query points from links
                    u_linked_du = patch.deriv_links_(:, 4*jj-3); %query points derivatives from links
                    u_linked_dv = patch.deriv_links_(:, 4*jj-2); %query points derivatives from links
                    v_linked_du = patch.deriv_links_(:, 4*jj-1); %query points derivatives from links
                    v_linked_dv = patch.deriv_links_(:, 4*jj); %query points derivatives from links
                    %vq = patch.deriv_links(:, 2*jj); %query points derivatives from links
                    f_du_linked_interp = interpPatch3(f_du_linked_patch, xq, linked_patch);
                    f_dv_linked_interp = interpPatch3(f_dv_linked_patch, xq, linked_patch);
                    %[u_linked_du, u_linked_dv] = patch.grad_FDM(uq);
                    %[v_linked_du, v_linked_dv] = patch.grad_FDM(vq);
                    %isnan(f_du_linked_interp)
                    %isnan(f_dv_linked_interp)
                    if jj==5
                        test_du = (f_du_linked_interp.*u_linked_du + f_dv_linked_interp.*v_linked_du);
                        test_dv = (f_du_linked_interp.*u_linked_dv + f_dv_linked_interp.*v_linked_dv);
                        test_du(12)
                        test_dv(12)
                        %f_du_linked_interp(12)
                        %f_dv_linked_interp(12)
                    end
                    ppou = patch.pou_links_(:,jj);
                    f_du(:, ii) = f_du(:, ii) + (f_du_linked_interp.*u_linked_du + f_dv_linked_interp.*v_linked_du).*patch.pou_links_(:,jj);
                    f_dv(:, ii) = f_dv(:, ii) + (f_du_linked_interp.*u_linked_dv + f_dv_linked_interp.*v_linked_dv).*patch.pou_links_(:,jj);  
                    disp('hello');
               end
            end
            

          end
          
      end

      function [f_du, f_dv] = patchwiseDerivFDM(obj, f)
          %Doesn't work for boundary of patch 
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          m = obj.patches(1).Nu; n = obj.patches(1).Nv;
          f_du = zeros((m)*(n), length(obj.patches));
          f_dv = zeros((m)*(n), length(obj.patches));
    
          for ii=1:length(obj.patches)
              patch = obj.patches(ii);
              f_patch = f(:, patch.numPatch);
              %fphi_patch = f_patch.*patch.pou;
%               if ii==pn
%                  stem3(reshape(f_patch,[m+2,n+2])); 
%               end
              [f_du_patch, f_dv_patch] = patch.grad_FDM(f_patch);

       
              f_du(:, patch.numPatch) = f_du_patch;
              f_dv(:, patch.numPatch) = f_dv_patch;
 

              
          end         
      end

      function f_blend = blendSurfaceFunction(obj, f)% tval1) %tval1 is true val on patch 1 
      % f is numNodes_*numPatches column vector function on all nodes of surface    
          numNodes = obj.patches(1).numNodes;
          f = reshape(f, [numNodes, obj.numPatches]);
          f_blend = zeros(numNodes, obj.numPatches);
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
              f_patch = f(:,ii);
              f_blend(:, ii) = f_patch.*patch.pou;
              
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);

                    f_linked_patch = f(:, jj);
                    xq = patch.cl_links(:, 2*jj-1:2*jj); %query points from links
                    
                    % ----check without interpolation error using exact values ----
                    %lap_phi = @(x,y,z) 8-3*(4*x.^2+4*y.^2);
                    %f_linked_interp = lap_phi(patch.x, patch.y, patch.z);
                    %  ---------------------------------------
                    % ----check without interpolation error using exact values on patch1 given in tval1 ----
                    %if ii==1
                        
                    %    f_linked_interp = tval1;
                    
                    %else 
                    %    f_linked_interp = interpPatch4(f_linked_patch, xq, linked_patch);
                    %end
                    %  ---------------------------------------                    
                    f_linked_interp = interpPatch4(f_linked_patch, xq, linked_patch);

                    
                    f_blend(:,ii) = f_blend(:,ii) + f_linked_interp.*patch.pou_links(:,jj);
                     
                    
               end
            end


          end
          
          f_blend = f_blend(:);
      end

      function normals_blend = getNormals(obj)
          %Function to get normals to the surface
            patch = obj.patches;
            x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
            y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
            z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];

            [x_du_app, x_dv_app] = obj.patchwiseDerivFDM(x);
            [y_du_app, y_dv_app] = obj.patchwiseDerivFDM(y);
            [z_du_app, z_dv_app] = obj.patchwiseDerivFDM(z);


            x_du = x_du_app(:);
            y_du = y_du_app(:);
            z_du = z_du_app(:);

            x_dv = x_dv_app(:);
            y_dv = y_dv_app(:);
            z_dv = z_dv_app(:);

            r_du = [x_du y_du z_du];
            r_dv = [x_dv y_dv z_dv];

            temp_cross = cross(r_du, r_dv, 2);
            vecn = vecnorm(temp_cross, 2, 2).^(-1);
            D = spdiags(vecn(:),0,length(vecn),length(vecn));
            normals_fdm_new = D*temp_cross;




            normals_fdm_new_blend_x = obj.blendSurfaceFunction(normals_fdm_new(:,1));
            normals_fdm_new_blend_y = obj.blendSurfaceFunction(normals_fdm_new(:,2));
            normals_fdm_new_blend_z = obj.blendSurfaceFunction(normals_fdm_new(:,3));
            normals_blend = [normals_fdm_new_blend_x normals_fdm_new_blend_y normals_fdm_new_blend_z ];

            %calculate normals using fdm          
      end
      
      
      %function to get mean  curvature and gaussian curvature of surface
      function [H_blend, K_blend] = getCurvature(obj)
            patch = obj.patches;
            x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
            y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
            z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];
    
    
            [x_du_app, x_dv_app] = obj.patchwiseDerivFDM(x);
            [y_du_app, y_dv_app] = obj.patchwiseDerivFDM(y);
            [z_du_app, z_dv_app] = obj.patchwiseDerivFDM(z);

            x_du = x_du_app(:);
            y_du = y_du_app(:);
            z_du = z_du_app(:);

            x_dv = x_dv_app(:);
            y_dv = y_dv_app(:);
            z_dv = z_dv_app(:);

            E = x_du.*x_du + y_du.*y_du + z_du.*z_du;
            F = x_du.*x_dv + y_du.*y_dv + z_du.*z_dv;
            G = x_dv.*x_dv + y_dv.*y_dv + z_dv.*z_dv;
            W = (E.*G - F.^2).^(0.5);


            r_du = [x_du y_du z_du];
            r_dv = [x_dv y_dv z_dv];

            temp_cross = cross(r_du, r_dv, 2);
            vecn = vecnorm(temp_cross, 2, 2).^(-1);
            D = spdiags(vecn(:),0,length(vecn),length(vecn));
            normals_fdm_new = D*temp_cross; %normals calculated

            normals_fdm_new_blend_x = obj.blendSurfaceFunction(normals_fdm_new(:,1));
            normals_fdm_new_blend_y = obj.blendSurfaceFunction(normals_fdm_new(:,2));
            normals_fdm_new_blend_z = obj.blendSurfaceFunction(normals_fdm_new(:,3));
            normals_fdm_new_blend = [normals_fdm_new_blend_x normals_fdm_new_blend_y normals_fdm_new_blend_z ];    

            [x_duu_app, x_duv_app] = obj.patchwiseDerivFDM(x_du_app);
            [x_dvu_app, x_dvv_app] = obj.patchwiseDerivFDM(x_dv_app);
            [y_duu_app, y_duv_app] = obj.patchwiseDerivFDM(y_du_app);
            [y_dvu_app, y_dvv_app] = obj.patchwiseDerivFDM(y_dv_app);    
            [z_duu_app, z_duv_app] = obj.patchwiseDerivFDM(z_du_app);
            [z_dvu_app, z_dvv_app] = obj.patchwiseDerivFDM(z_dv_app);

            

            x_duu = x_duu_app(:);
            x_duv = x_duv_app(:);
            x_dvv = x_dvv_app(:);
            y_duu = y_duu_app(:);
            y_duv = y_duv_app(:);
            y_dvv = y_dvv_app(:);
            z_duu = z_duu_app(:);
            z_duv = z_duv_app(:);
            z_dvv = z_dvv_app(:);

            nx = normals_fdm_new_blend(:,1);
            ny = normals_fdm_new_blend(:,2);
            nz = normals_fdm_new_blend(:,3);


            L = x_duu.*nx + y_duu.*ny + z_duu.*nz;
            M = x_duv.*nx + y_duv.*ny + z_duv.*nz;
            N = x_dvv.*nx + y_dvv.*ny + z_dvv.*nz;



            H = -0.5*(E.*N - 2*F.*M + G.*L)./ W.^2;
            H_blend = obj.blendSurfaceFunction(H); 
            
            K = (L.*N - M.^2)./ W.^2;
            K_blend = obj.blendSurfaceFunction(K); 
                
          
          
      end
      
      
      function surf_lap_blend = getLaplacian(obj, f)
      %Function to get surface laplacian of f 
            f = reshape(f, [obj.patches(1).numNodes, obj.patches(1).numPatches]);
            patch = obj.patches;
            x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
            y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
            z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];


            [x_du_app, x_dv_app] = obj.patchwiseDerivFDM(x);
            [y_du_app, y_dv_app] = obj.patchwiseDerivFDM(y);
            [z_du_app, z_dv_app] = obj.patchwiseDerivFDM(z);

            x_du = x_du_app(:);
            y_du = y_du_app(:);
            z_du = z_du_app(:);

            x_dv = x_dv_app(:);
            y_dv = y_dv_app(:);
            z_dv = z_dv_app(:);



            [f_du_app, f_dv_app] = obj.patchwiseDerivFDM(f);

            f_du = f_du_app(:);
            f_dv = f_dv_app(:);


            E = x_du.*x_du + y_du.*y_du + z_du.*z_du;
            F = x_du.*x_dv + y_du.*y_dv + z_du.*z_dv;
            G = x_dv.*x_dv + y_dv.*y_dv + z_dv.*z_dv;
            W = (E.*G - F.^2).^(0.5);

            % calculate ((E*f_dv - F*f_du)/W) and take v derivative 
            t1 = (E.*f_dv - F.*f_du)./W;
            t1_app = reshape(t1, [size(x,1),size(x,2)]);

            [t1_du_app, t1_dv_app] = obj.patchwiseDerivFDM(t1_app);
            t1_du = t1_du_app(:);
            t1_dv = t1_dv_app(:);


            % calculate ((G*f_du - F*f_dv)/W) and take u derivative

            t2 = (G.*f_du - F.*f_dv)./W;
            t2_app = reshape(t2, [size(x,1),size(x,2)]);

            [t2_du_app, t2_dv_app] = obj.patchwiseDerivFDM(t2_app);
            t2_du = t2_du_app(:);
            t2_dv = t2_dv_app(:);

            %calculate surf_lap

            surf_lap = (t1_dv + t2_du)./W;

            surf_lap_blend = obj.blendSurfaceFunction(surf_lap);
%             true_lap = true_lap(:);
% 
%             error_lap = max(abs(surf_lap-true_lap))/max(abs(true_lap))
%             error_lap_blend = max(abs(surf_lap_blend-true_lap))/max(abs(true_lap))

          
          
      end

      function surf_grad_blend = getSurfaceGradient(obj, f)
      %Function to get surface gradient of scalar f 
            f = reshape(f, [obj.patches(1).numNodes, obj.patches(1).numPatches]);
            patch = obj.patches;
            x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
            y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
            z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];


            [x_du_app, x_dv_app] = obj.patchwiseDerivFDM(x);
            [y_du_app, y_dv_app] = obj.patchwiseDerivFDM(y);
            [z_du_app, z_dv_app] = obj.patchwiseDerivFDM(z);

            x_du = x_du_app(:);
            y_du = y_du_app(:);
            z_du = z_du_app(:);

            x_dv = x_dv_app(:);
            y_dv = y_dv_app(:);
            z_dv = z_dv_app(:);



            [f_du_app, f_dv_app] = obj.patchwiseDerivFDM(f);

            f_du = f_du_app(:);
            f_dv = f_dv_app(:);


            E = x_du.*x_du + y_du.*y_du + z_du.*z_du;
            F = x_du.*x_dv + y_du.*y_dv + z_du.*z_dv;
            G = x_dv.*x_dv + y_dv.*y_dv + z_dv.*z_dv;
            W = (E.*G - F.^2).^(0.5);

            % calculate ((E*f_dv - F*f_du)/W) and take v derivative 
            t1_x = ((G.*x_du - F.*x_dv)./W.^2).*f_du;
            t1_y = ((G.*y_du - F.*y_dv)./W.^2).*f_du;
            t1_z = ((G.*z_du - F.*z_dv)./W.^2).*f_du;
            t2_x = ((E.*x_dv - F.*x_du)./W.^2).*f_dv;
            t2_y = ((E.*y_dv - F.*y_du)./W.^2).*f_dv;
            t2_z = ((E.*z_dv - F.*z_du)./W.^2).*f_dv;
                  

            %calculate surf_grad

            surf_grad = [t1_x+t2_x t1_y+t2_y t1_z+t2_z];

            surf_grad_blend_x = obj.blendSurfaceFunction(surf_grad(:,1));
            surf_grad_blend_y = obj.blendSurfaceFunction(surf_grad(:,2));
            surf_grad_blend_z = obj.blendSurfaceFunction(surf_grad(:,3));
            surf_grad_blend = [surf_grad_blend_x surf_grad_blend_y surf_grad_blend_z];
          
      end      

      function fb = getBendingForce(obj)
         %get bending force 
         
         [H, K] = obj.getCurvature();
         n = obj.getNormals();
         lapH = obj.getLaplacian(H);
         temp = obj.kb*(lapH + 2*H.*(H.^2-K));
         temp = repmat(temp, [1,3]);
         fb = temp.*n;
         
          
          
      end

      
      function f_sigma = getTensileForce(obj, sigma)
          %get tensile force
            patch = obj.patches;
            x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
            y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
            z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];
          
            lap_x = obj.getLaplacian(x(:));
            lap_y = obj.getLaplacian(y(:));
            lap_z = obj.getLaplacian(z(:));
            lap = [lap_x lap_y lap_z];
            sigma_ = repmat(sigma, [1,3]);
            
            sigma_lap = sigma_.*lap;
            
            grad_sigma = obj.getSurfaceGradient(sigma);
            
            f_sigma = sigma_lap + grad_sigma;
            
%             [H,K] = obj.getCurvature();
%             n = obj.getNormals();
%             H_ = repmat(H,[1,3]);
%             f_sigma = 2*sigma_.*H_.*n + grad_sigma;
            
            
                
      end
      
      function val = integrateOverSurface(obj, f)
          %function to integrate f over surface
          %f is column vector, 6*m*n x 1 size. 
          val = 0;
          f_mat = reshape(f, [obj.patches(1).Nu*obj.patches(1).Nv, length(obj.patches)]);
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              val = val + integratePatch(patch, f_mat(:,ii));
               
          end   
          
          
          
      end
      
      function area = getArea(obj)
          f = ones(obj.patches(1).Nu*obj.patches(1).Nv*length(obj.patches), 1);
          area = integrateOverSurface(obj, f);
          
      end

      function vol = getVolume(obj)
          patch = obj.patches;
          x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
          x = x(:);
          f = zeros(size(x, 1), 3);
          f(:,1) = x;
          normals = obj.getNormals();
          f = sum(f.*normals, 2);
          vol = integrateOverSurface(obj, f);
          
      end      

      function I = getInertia(obj)
          I = zeros(3,3);
          patch = obj.patches;
          
          x = [patch(1).r(:,1) patch(2).r(:,1) patch(3).r(:,1) patch(4).r(:,1) patch(5).r(:,1) patch(6).r(:,1)];
          y = [patch(1).r(:,2) patch(2).r(:,2) patch(3).r(:,2) patch(4).r(:,2) patch(5).r(:,2) patch(6).r(:,2)];
          z = [patch(1).r(:,3) patch(2).r(:,3) patch(3).r(:,3) patch(4).r(:,3) patch(5).r(:,3) patch(6).r(:,3)];
                    
          x = x(:); y = y(:); z = z(:);
          f_xx = zeros(size(x, 1), 3);
          f_xy = zeros(size(x, 1), 3);
          f_xz = zeros(size(x, 1), 3);
          f_yy = zeros(size(x, 1), 3);
          f_yz = zeros(size(x, 1), 3);
          f_zz = zeros(size(x, 1), 3);
          
          f_xx(:,1) = x.*y.^2; f_xx(:,2) = y.*z.^2;
          f_xy(:,3) = -x.*y.*z; 
          f_xz(:,2) = -x.*y.*z;
          f_yy(:,1) = x.*z.^2; f_yy(:,2) = y.*x.^2;
          f_yz(:,1) = -x.*y.*z;
          f_zz(:,2) = y.*x.^2; f_zz(:,3) = z.*y.^2;
          
          normals = obj.getNormals();
          f_xx = sum(f_xx.*normals, 2);
          f_xy = sum(f_xy.*normals, 2);
          f_xz = sum(f_xz.*normals, 2);
          f_yy = sum(f_yy.*normals, 2);
          f_yz = sum(f_yz.*normals, 2);
          f_zz = sum(f_zz.*normals, 2);
          I(1,1) = integrateOverSurface(obj, f_xx);
          I(1,2) = integrateOverSurface(obj, f_xy);
          I(1,3) = integrateOverSurface(obj, f_xz);
          I(2,1) = I(1,2);
          I(2,2) = integrateOverSurface(obj, f_yy);
          I(2,3) = integrateOverSurface(obj, f_yz);
          I(3,1) = I(1,3);
          I(3,2) = I(2,3);
          I(3,3) = integrateOverSurface(obj, f_zz);
          
      end            

      function [f_du, f_dv] = unique_uv_deriv(obj, f)
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          f_du = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dv = zeros(obj.patches(1).numNodes, length(obj.patches));
          %f_du_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
          %f_dv_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
         
          [f_du_patches, f_dv_patches] = obj.patchwiseDerivFDM(f);
          %f_du = f_du_patches; f_dv = f_dv_patches;

          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
              f_du(:, ii)   = f_du_patches(:, ii).*patch.pou;
              f_dv(:, ii)   = f_dv_patches(:, ii).*patch.pou;
              
          end
          
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    f_du_linked_patch = f_du_patches(:, jj);
                    f_dv_linked_patch = f_dv_patches(:, jj);
                    xq = patch.links(:, 2*jj-1:2*jj); %query points from links
                    u_linked_du = patch.deriv_links(:, 4*jj-3); %query points derivatives from links
                    u_linked_dv = patch.deriv_links(:, 4*jj-2); %query points derivatives from links
                    v_linked_du = patch.deriv_links(:, 4*jj-1); %query points derivatives from links
                    v_linked_dv = patch.deriv_links(:, 4*jj); %query points derivatives from links
                    %vq = patch.deriv_links(:, 2*jj); %query points derivatives from links
                    f_du_linked_interp = interpPatch4(f_du_linked_patch, xq, linked_patch);
                    f_dv_linked_interp = interpPatch4(f_dv_linked_patch, xq, linked_patch);
                    %[u_linked_du, u_linked_dv] = patch.grad_FDM(uq);
                    %[v_linked_du, v_linked_dv] = patch.grad_FDM(vq);
                    f_du(:, ii) = f_du(:, ii) + (f_du_linked_interp.*u_linked_du + f_dv_linked_interp.*v_linked_du).*patch.pou_links(:, jj);
                    f_dv(:, ii) = f_dv(:, ii) + (f_du_linked_interp.*u_linked_dv + f_dv_linked_interp.*v_linked_dv).*patch.pou_links(:, jj);  
                    
               end
            end
            

          end
      end      

      
      function [f_duu, f_dvv] = unique_uv_second_deriv(obj, f)
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          f_duu = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_duv = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dvu = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dvv = zeros(obj.patches(1).numNodes, length(obj.patches));
          %f_du_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
          %f_dv_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
         
          [f_du_patches, f_dv_patches] = obj.patchwiseDerivFDM(f);
          [f_duu_patches, f_duv_patches] = obj.patchwiseDerivFDM(f_du_patches);
          [f_dvu_patches, f_dvv_patches] = obj.patchwiseDerivFDM(f_dv_patches);
          %f_du = f_du_patches; f_dv = f_dv_patches;

          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
              f_duu(:, ii)   = f_duu_patches(:, ii).*patch.pou;
              f_duv(:, ii)   = f_duv_patches(:, ii).*patch.pou;
              f_dvu(:, ii)   = f_dvu_patches(:, ii).*patch.pou;
              f_dvv(:, ii)   = f_dvv_patches(:, ii).*patch.pou;
                       
          end
          
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    f_du_linked_patch = f_du_patches(:, jj);
                    f_dv_linked_patch = f_dv_patches(:, jj);
                    f_duu_linked_patch = f_duu_patches(:, jj);
                    f_duv_linked_patch = f_duv_patches(:, jj);
                    f_dvu_linked_patch = f_dvu_patches(:, jj);
                    f_dvv_linked_patch = f_dvv_patches(:, jj);                    
                    xq = patch.links(:, 2*jj-1:2*jj); %query points from links
                    u_linked_du = patch.deriv_links(:, 4*jj-3); %query points derivatives from links
                    u_linked_dv = patch.deriv_links(:, 4*jj-2); %query points derivatives from links
                    v_linked_du = patch.deriv_links(:, 4*jj-1); %query points derivatives from links
                    v_linked_dv = patch.deriv_links(:, 4*jj); %query points derivatives from links
                    
                    u_linked_duu = patch.secondderiv_links(:, 8*jj-7); %query points second derivatives from links
                    u_linked_duv = patch.secondderiv_links(:, 8*jj-6); %query points second derivatives from links
                    u_linked_dvu = patch.secondderiv_links(:, 8*jj-5); %query points second derivatives from links
                    u_linked_dvv = patch.secondderiv_links(:, 8*jj-4); %query points second derivatives from links
                    v_linked_duu = patch.secondderiv_links(:, 8*jj-3); %query points second derivatives from links
                    v_linked_duv = patch.secondderiv_links(:, 8*jj-2); %query points second derivatives from links
                    v_linked_dvu = patch.secondderiv_links(:, 8*jj-1); %query points second derivatives from links
                    v_linked_dvv = patch.secondderiv_links(:, 8*jj); %query points second derivatives from links
                                  
                    %vq = patch.deriv_links(:, 2*jj); %query points derivatives from links
                    f_du_linked_interp = interpPatch4(f_du_linked_patch, xq, linked_patch);
                    f_dv_linked_interp = interpPatch4(f_dv_linked_patch, xq, linked_patch);                    
                    f_duu_linked_interp = interpPatch4(f_duu_linked_patch, xq, linked_patch);
                    f_duv_linked_interp = interpPatch4(f_duv_linked_patch, xq, linked_patch);
                    f_dvu_linked_interp = interpPatch4(f_dvu_linked_patch, xq, linked_patch);
                    f_dvv_linked_interp = interpPatch4(f_dvv_linked_patch, xq, linked_patch);                    
                    %[u_linked_du, u_linked_dv] = patch.grad_FDM(uq);
                    %[v_linked_du, v_linked_dv] = patch.grad_FDM(vq);
                    f_duu(:, ii) = f_duu(:, ii) +  (f_duu_linked_interp.*u_linked_du + f_dvu_linked_interp.*v_linked_du).*u_linked_du.*patch.pou_links(:, jj) + (f_duv_linked_interp.*u_linked_du + f_dvv_linked_interp.*v_linked_du).*v_linked_du.*patch.pou_links(:, jj) + (f_du_linked_interp.*u_linked_duu + f_dv_linked_interp.*v_linked_duu).*patch.pou_links(:, jj);
                    f_dvv(:, ii) = f_dvv(:, ii) +  (f_duu_linked_interp.*u_linked_dv + f_dvu_linked_interp.*v_linked_dv).*u_linked_dv.*patch.pou_links(:, jj) + (f_duv_linked_interp.*u_linked_dv + f_dvv_linked_interp.*v_linked_dv).*v_linked_dv.*patch.pou_links(:, jj) + (f_du_linked_interp.*u_linked_dvv + f_dv_linked_interp.*v_linked_dvv).*patch.pou_links(:, jj);
                       
                    
                    
               end
            end
            

          end
      end       
      
      function [f_du_chain, f_dv_chain] = chain_rule(obj, f_du, f_dv)
          %Function to calculate derivative of scalar f on surface S
          % f is to be arranged in order of patch nodes, each patch stacked
          % like a column, # of cols = # of patches
          f_du_chain = zeros(obj.patches(1).numNodes, length(obj.patches));
          f_dv_chain = zeros(obj.patches(1).numNodes, length(obj.patches));
          %f_du_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
          %f_dv_patches = zeros(obj.patches(1).numNodes, length(obj.patches));
         
          %[f_du_patches, f_dv_patches] = obj.patchwiseDerivFDM(f);
          

          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
              f_du_chain(:, ii)   = f_du(:, ii).*patch.pou;
              f_dv_chain(:, ii)   = f_dv(:, ii).*patch.pou;
              
          end
          
          for ii = 1:obj.numPatches
              patch = obj.patches(ii);
            for jj=1:obj.numPatches
                  
               if jj ~= ii
                    %f_linked_patch = f(:, jj);
                    linked_patch = obj.patches(jj);
                    f_du_linked_patch = f_du(:, jj);
                    f_dv_linked_patch = f_dv(:, jj);
                    xq = patch.links(:, 2*jj-1:2*jj); %query points from links
                    u_linked_du = patch.deriv_links(:, 4*jj-3); %query points derivatives from links
                    u_linked_dv = patch.deriv_links(:, 4*jj-2); %query points derivatives from links
                    v_linked_du = patch.deriv_links(:, 4*jj-1); %query points derivatives from links
                    v_linked_dv = patch.deriv_links(:, 4*jj); %query points derivatives from links
                    %vq = patch.deriv_links(:, 2*jj); %query points derivatives from links
                    f_du_linked_interp = interpPatch4(f_du_linked_patch, xq, linked_patch);
                    f_dv_linked_interp = interpPatch4(f_dv_linked_patch, xq, linked_patch);
                    %[u_linked_du, u_linked_dv] = patch.grad_FDM(uq);
                    %[v_linked_du, v_linked_dv] = patch.grad_FDM(vq);
                    f_du_chain(:, ii) = f_du_chain(:, ii) + (f_du_linked_interp.*u_linked_du + f_dv_linked_interp.*v_linked_du).*patch.pou_links(:, jj);
                    f_dv_chain(:, ii) = f_dv_chain(:, ii) + (f_du_linked_interp.*u_linked_dv + f_dv_linked_interp.*v_linked_dv).*patch.pou_links(:, jj);  
                    
               end
            end
            

          end
      end      
      
      % function p2p_deriv: nodes of p1 transformed to u-v coords in p2 and then calculate derivative
      % of f (f is provided on p2) in p2 and then interpolate to tranformed
      % u-v in p2. 
       function [f_du_interp, f_dv_interp] = p2p_deriv(obj, f, p1, p2)
          %Function to calculate derivative of scalar f on patch p2 at
          %linked nodes of patch p1 
          
         
          %f = f.*(p2.pou);
          [f_du, f_dv] = p2.grad_patch(f);
          jj = p2.numPatch;
          xq = p1.links(:, 2*jj-1:2*jj); %query points from links
          f_du_interp = interpPatch2(f_du, xq, p2);
          f_dv_interp = interpPatch2(f_dv, xq, p2);
          
         
          
  
       end
      
       
       function surf_lap_blend = getSphereLaplacian(obj, f)
      %Function to get sphere surface laplacian of f 
            f = reshape(f, [obj.patches(1).numNodes, obj.patches(1).numPatches]);
            patch = obj.patches;
            x = [patch(1).sphcart(:,1) patch(2).sphcart(:,1) patch(3).sphcart(:,1) patch(4).sphcart(:,1) patch(5).sphcart(:,1) patch(6).sphcart(:,1)];
            y = [patch(1).sphcart(:,2) patch(2).sphcart(:,2) patch(3).sphcart(:,2) patch(4).sphcart(:,2) patch(5).sphcart(:,2) patch(6).sphcart(:,2)];
            z = [patch(1).sphcart(:,3) patch(2).sphcart(:,3) patch(3).sphcart(:,3) patch(4).sphcart(:,3) patch(5).sphcart(:,3) patch(6).sphcart(:,3)];


            [x_du_app, x_dv_app] = obj.patchwiseDerivFDM(x);
            [y_du_app, y_dv_app] = obj.patchwiseDerivFDM(y);
            [z_du_app, z_dv_app] = obj.patchwiseDerivFDM(z);

            x_du = x_du_app(:);
            y_du = y_du_app(:);
            z_du = z_du_app(:);

            x_dv = x_dv_app(:);
            y_dv = y_dv_app(:);
            z_dv = z_dv_app(:);



            [f_du_app, f_dv_app] = obj.patchwiseDerivFDM(f);

            f_du = f_du_app(:);
            f_dv = f_dv_app(:);


            E = x_du.*x_du + y_du.*y_du + z_du.*z_du;
            F = x_du.*x_dv + y_du.*y_dv + z_du.*z_dv;
            G = x_dv.*x_dv + y_dv.*y_dv + z_dv.*z_dv;
            W = (E.*G - F.^2).^(0.5);

            % calculate ((E*f_dv - F*f_du)/W) and take v derivative 
            t1 = (E.*f_dv - F.*f_du)./W;
            t1_app = reshape(t1, [size(x,1),size(x,2)]);

            [t1_du_app, t1_dv_app] = obj.patchwiseDerivFDM(t1_app);
            t1_du = t1_du_app(:);
            t1_dv = t1_dv_app(:);


            % calculate ((G*f_du - F*f_dv)/W) and take u derivative

            t2 = (G.*f_du - F.*f_dv)./W;
            t2_app = reshape(t2, [size(x,1),size(x,2)]);

            [t2_du_app, t2_dv_app] = obj.patchwiseDerivFDM(t2_app);
            t2_du = t2_du_app(:);
            t2_dv = t2_dv_app(:);

            %calculate surf_lap

            surf_lap = (t1_dv + t2_du)./W;

            surf_lap_blend = obj.blendSurfaceFunction(surf_lap);
%             true_lap = true_lap(:);
% 
%             error_lap = max(abs(surf_lap-true_lap))/max(abs(true_lap))
%             error_lap_blend = max(abs(surf_lap_blend-true_lap))/max(abs(true_lap))

          
          
       end
      
      function surf_grad_blend = getSphereGradient(obj, f)
      %Function to get surface gradient of scalar f 
            f = reshape(f, [obj.patches(1).numNodes, obj.patches(1).numPatches]);
            patch = obj.patches;
            x = [patch(1).sphcart(:,1) patch(2).sphcart(:,1) patch(3).sphcart(:,1) patch(4).sphcart(:,1) patch(5).sphcart(:,1) patch(6).sphcart(:,1)];
            y = [patch(1).sphcart(:,2) patch(2).sphcart(:,2) patch(3).sphcart(:,2) patch(4).sphcart(:,2) patch(5).sphcart(:,2) patch(6).sphcart(:,2)];
            z = [patch(1).sphcart(:,3) patch(2).sphcart(:,3) patch(3).sphcart(:,3) patch(4).sphcart(:,3) patch(5).sphcart(:,3) patch(6).sphcart(:,3)];


            [x_du_app, x_dv_app] = obj.patchwiseDerivFDM(x);
            [y_du_app, y_dv_app] = obj.patchwiseDerivFDM(y);
            [z_du_app, z_dv_app] = obj.patchwiseDerivFDM(z);

            x_du = x_du_app(:);
            y_du = y_du_app(:);
            z_du = z_du_app(:);

            x_dv = x_dv_app(:);
            y_dv = y_dv_app(:);
            z_dv = z_dv_app(:);



            [f_du_app, f_dv_app] = obj.patchwiseDerivFDM(f);

            f_du = f_du_app(:);
            f_dv = f_dv_app(:);


            E = x_du.*x_du + y_du.*y_du + z_du.*z_du;
            F = x_du.*x_dv + y_du.*y_dv + z_du.*z_dv;
            G = x_dv.*x_dv + y_dv.*y_dv + z_dv.*z_dv;
            W = (E.*G - F.^2).^(0.5);

            % calculate ((E*f_dv - F*f_du)/W) and take v derivative 
            t1_x = ((G.*x_du - F.*x_dv)./W.^2).*f_du;
            t1_y = ((G.*y_du - F.*y_dv)./W.^2).*f_du;
            t1_z = ((G.*z_du - F.*z_dv)./W.^2).*f_du;
            t2_x = ((E.*x_dv - F.*x_du)./W.^2).*f_dv;
            t2_y = ((E.*y_dv - F.*y_du)./W.^2).*f_dv;
            t2_z = ((E.*z_dv - F.*z_du)./W.^2).*f_dv;
                  

            %calculate surf_grad

            surf_grad = [t1_x+t2_x t1_y+t2_y t1_z+t2_z];

            surf_grad_blend_x = obj.blendSurfaceFunction(surf_grad(:,1));
            surf_grad_blend_y = obj.blendSurfaceFunction(surf_grad(:,2));
            surf_grad_blend_z = obj.blendSurfaceFunction(surf_grad(:,3));
            surf_grad_blend = [surf_grad_blend_x surf_grad_blend_y surf_grad_blend_z];
          
      end      
   end
      
 
end


