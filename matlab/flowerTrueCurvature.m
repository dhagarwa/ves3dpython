function [H_true,K_true, n_true, lapH_true] = flowerTrueCurvature(m)
%calculate true mean curvature of a flower
%caluclate true curvature
    clc
    n = m;
    syms th ph
    %th = u, ph = v
    rho_ =  1 + exp(-3*(1/4)*sqrt(105/(2*pi))*sin(th)^2*cos(th)*cos(2*ph));
    x_ = rho_*sin(th)*cos(ph); y_ = rho_*sin(th)*sin(ph); z_ = rho_*cos(th); 

    X  = [x_;y_;z_]; 
    Xu = (diff(X,th));
    Xv = (diff(X,ph));

    Xuu = (diff(X,th,2));
    Xuv = (diff(Xv,th)); 
    Xvv = (diff(X,ph,2)); 
    disp('....');

    [the,phi]=ndgrid(pi*(1:m)/(m+1), (pi )*(1:n)/(n+1));
    the = the(:); phi = phi(:);
    the = the(1:2); phi = phi(1:2);
    %size(phi)
    %Xu = double(subs(Xu, {'th' 'ph'}, {the phi}));
    %Xv = double(subs(Xv, {'th' 'ph'}, {the phi}));
    %Xuu = double(subs(Xuu, {'th' 'ph'}, {the phi}));
    %Xuv = double(subs(Xuv, {'th' 'ph'}, {the phi}));
    %Xvv = double(subs(Xvv, {'th' 'ph'}, {the phi}));
    
    E_ = (Xu.'*Xu);
    F_ = (Xu.'*Xv);
    G_ = (Xv.'*Xv);
    W_ = sqrt(E_*G_ - F_^2);
    disp('......');

    nor = cross(Xu,Xv)/W_;
    disp('........');

    L_ = (Xuu.'*nor);
    M_ = (Xuv.'*nor);
    N_ = (Xvv.'*nor);
    disp('..........');

    H_ = -(((E_*N_ - 2*F_*M_ + G_*L_))/W_^2/2);
    K_ = (((L_*N_ - M_^2)/W_^2));
    % laplacian of H
    H_u = diff(H_,th);
    H_v = diff(H_,ph);
    t1 = (E_*H_v-F_*H_u)/W_; t1_v = diff(t1,ph);
    t2 = (G_*H_u-F_*H_v)/W_; t2_u = diff(t2,th);
    lapH = (t1_v+t2_u)/W_ ; 
    
    
    H_true = double(subs(H_, {'th' 'ph'}, {the phi}));
    K_true = double(subs(K_, {'th' 'ph'}, {the phi}));
    n_true = double(subs(nor, {'th' 'ph'}, {the phi}));
    %size(n_true)
    lapH_true = double(subs(lapH, {'th' 'ph'}, {the phi}));
    
    
    
    
    
    
end