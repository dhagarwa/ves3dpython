classdef Surface < handle
   properties
        patches
        link %1 linked , 0 not linked
   end
   
   methods
      function obj = Surface(patches)
            obj.patches = patches;
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
         obj.patches(1).links
          val = 1;
      end
      
      
   end
      
 
end


