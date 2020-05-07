classdef Surface < handle
   properties
        patches
        link %1 linked , 0 not linked
        numPatches
   end
   
   methods
      function obj = Surface(patches)
            obj.patches = patches;
            obj.numPatches = length(patches);
            obj.link = obj.linkPatches();
            
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


