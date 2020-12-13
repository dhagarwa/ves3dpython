function [H_true,n_true] = flowerTrueCurvature(m)
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
    H_true = double(subs(H_, {'th' 'ph'}, {the phi}));
    n_true = double(subs(nor, {'th' 'ph'}, {the phi}));

    %K = simplify(((L*N - M^2)/W^2));
    
    
    
    
    
end