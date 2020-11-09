function [r1, r2] = testlinks(S)
    patches = S.patches;
    err = [];
    err2 = [];
    data = [];
    for n1 = 1:6
        for n2 = 1:6
    p1 = patches(n1);
    p2 = patches(n2);

    for ii=1:size(p1.r, 1)
       
        r1 = p1.r(ii, :);
        node2 = patchParameterise(r1, p1, p2);
        u1 = p1.u(ii);
        v1 = p1.v(ii);
        u2 = node2(1);
        v2 = node2(2);
        r2 = sph2cartPatch(u2, v2, p2);
        %d_r = [r1; r2 ]
        %d_u = [u1 v1; u2 v2]
        err = [err; norm(r2-r1)];
        %err2 = [err2; norm([u1 v1]-[u2 v2])];
    end
    %err
        end
    end
    max(err)
    %max(err2)
    

    

end