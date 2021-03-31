function [] = plotSurface(S)
    for ii=1:S.numPatches
        patch = S.patches(ii);
        tri = delaunay(patch.x, patch.y);
        %plot(x,y,'.')
%%
% How many triangles are there?
        [r,c] = size(tri);
        %disp(r)
%% Plot it with TRISURF
        h = trisurf(tri, patch.x, patch.y, patch.z);
        axis vis3d
        hold on
        axis([-1 1 -1 1 -1 1])
        
    end

end