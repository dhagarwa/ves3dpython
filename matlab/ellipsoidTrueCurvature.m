function [H_true,K_true, n_true, lapH_true] = ellipsoidTrueCurvature(m)
%calculate true mean curvature of an ellipsoid
%caluclate true curvature
    clc
    n = m;
    syms th ph
    %th = u, ph = v
    
    x_ = 1*sin(th)*cos(ph); y_ = 0.4*sin(th)*sin(ph); z_ = 1*cos(th); 

    X  = [x_;y_;z_]; 
    Xu = (diff(X,th));
    Xv = (diff(X,ph));

    Xuu = (diff(X,th,2));
    Xuv = (diff(Xv,th)); 
    Xvv = (diff(X,ph,2)); 
    disp('....');

    [the,phi]=ndgrid(pi*(1:m)/(m+1), (pi )*(1:n)/(n+1));
    the = the(:); phi = phi(:);
    the = the(1001:1050); phi = phi(1001:1050);
    
    %read theta phi from files 
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'theta_out32.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    the = cell2mat(C);
    
         
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'phi_out32.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    phi = cell2mat(C);
       
    
    
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'test_ellipsoid2_fb_src32.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    A = cell2mat(C);
    m0 = size(A, 1);
    fb = reshape(A, [m0/3,3]);
    %select few theta phi from all 
    the = the(1:10:1000); phi = phi(1:10:1000);
    fb = -fb(1:10:1000, :);
    %reading from files ends
    
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
    size(n_true);
    lapH_true = double(subs(lapH, {'th' 'ph'}, {the phi}));
    %max(-lapH_true)
    n_true = reshape(n_true, [100, 3]);
    fb_true = repmat(lapH_true + 2*H_true.*(H_true.^2 - K_true), [1, 3]).*n_true;
    
    error_fb_SPH = max(vecnorm(fb_true - fb, 2, 2))/max(vecnorm(fb_true,2,2))
    
    
    
    
    
end