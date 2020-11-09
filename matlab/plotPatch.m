function [] = plotPatch(S, numPatch)
        patch = S.patches(numPatch);
        tri = delaunay(patch.x, patch.y);
        %plot(x,y,'.')
%%
% How many triangles are there?
        [r,c] = size(tri);
        %disp(r)
%% Plot it with TRISURF
        h = trisurf(tri, patch.x, patch.y, patch.z);
        axis vis3d
        %hold on
        
    

end