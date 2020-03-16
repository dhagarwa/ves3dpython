function [r1, r2] = testlinks(S)
    patches = S.patches;
    patch = patches(1);

    u2 = patch.links(:, 5);
    v2 = patch.links(:, 6);
    r2 = sph2cartPatch(u2, v2, patches(3));
    r1 = patches.r;
    

end