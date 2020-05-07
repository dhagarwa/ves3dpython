function [r1, r2] = testlinks(S)
    patches = S.patches;
    patch = patches(5);

    u2 = patch.links(:, 7);
    v2 = patch.links(:, 8);
    r2 = sph2cartPatch(u2, v2, patches(4));
    r1 = patch.r;
    

end