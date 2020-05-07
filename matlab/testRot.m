function [] = testRot(S)
    p1 = S.patches(1);
    p2 = S.patches(3);
    n = 61;
    u1 = p1.u(n);
    v1 = p1.v(n);
    r1 = sph2cartPatch(u1, v1, p1)
    u2 = p1.links(n, 5);
    v2 = p1.links(n, 6);
    r2 = sph2cartPatch(u2, v2, p2)
    eps = 0.001;
    r1ue = sph2cartPatch(u1+eps, v1, p1);
    r2ue = sph2cartPatch(u2+eps, v2, p2);
    r1ve = sph2cartPatch(u1, v1+eps, p1);
    r2ve = sph2cartPatch(u2, v2+eps, p2);
    du1 = (r1ue - r1)/eps;
    dv1 = (r1ve - r1)/eps;
    du2 = (r2ue - r2)/eps;
    dv2 = (r2ve - r2)/eps;
    du1 = du1/norm(du1);
    dv1 = dv1/norm(dv1);
    du2 = du2/norm(du2);
    dv2 = dv2/norm(dv2);
    cross(du1, dv1)
    cross(du2, dv2)
    angle(du1, -du2)*2/pi
    angle(du1, -dv2)*2/pi
    R = [1 0 0; 0 0 1; 0 1 0];
    %norm(R*du2' - du1')
    %norm(R*dv2' - dv1')
    


end