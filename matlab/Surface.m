classdef Surface < handle
   properties
        patches
        link %1 linked , 0 not linked
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
          for ii=1:length(obj.patches)              
              patch = obj.patches(ii);
              f_ii = patch.shearForce(obj.Es, obj.Ed);
              patch.q_sl = f_ii;
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
                    x_linked_patch = linked_patch.pou;
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
               r0 = patch.r(jj, :);
               for kk=1:length(obj.patches)
                  p = obj.patches(kk);
                  if kk ~= patch.numPatch
                     %find out u-v coordinates in patch p 
                      node = patchParameterise(r0, patch, p);
                      patch.links(jj, 2*(kk-1) + 1:2*kk) = node;
                  end
               end
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
      
   end
      
 
end


