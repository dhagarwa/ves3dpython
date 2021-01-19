function [] = ellipsoidSLBendingTest4()
      clc

    m = 63;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'ellipsoid2_matlabtrg32.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    A = cell2mat(C);
    m0 = size(A, 1);
    all_trg = reshape(A, [3,m0/3]);
    all_trg = all_trg';

    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'test_ellipsoid2_SL_fb_matlabtrg32_p64.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    A = cell2mat(C);
    m0 = size(A, 1);
    all_pot = reshape(A, [3,m0/3]);
    all_pot = all_pot';
    
%     rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
%     relativefolder = 'testfiles';
%     fid = fopen(fullfile(rootdir, relativefolder, 'test_ellipsoid2_SL_fb_trg64.txt'), 'rt');
%     C = textscan(fid,'%f');
%     fclose(fid);
%     %celldisp(C);
%     A = cell2mat(C);
%     m0 = size(A, 1);
%     all_pot_trg = reshape(A, [3, m0/3]);
%     all_pot_trg = all_pot_trg';    
% 
%     error_src_trg = max(vecnorm(all_pot - all_pot_trg,2,2))/max(vecnorm(all_pot,2,2))
    
    patches = [];
    trg = [0.1690   -0.0903   -2.0063]; val = 0;
    for i=1:6
       patch =  standardEllipsoidPatch(m, n, i, R, 1, 0.4, 1);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = ones(size(patch.r,1), 3);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, 1, patch, patch.q_dl);
       
        
    end
%     mtrg = [];
%     for i=1:6
%        patch = patches(i);
%        mtrg = [mtrg; patch.r];
%     end
%     mtrg = mtrg';
%     mtrg = mtrg(:);
%     dlmwrite('flower_trg.txt',mtrg,'newline','pc');
%     
    S = Surface(patches, [1, 0, 0]);
    %DLvalue = DLSmooth(trg,1,S)
    
    %assign bending force as q_sl
    fb = S.getBendingForce();
    [H,K] = S.getCurvature();
    lapH = S.getLaplacian(H);
%     [H_true,K_true, n_true, lapH_true] = ellipsoidTrueCurvature(m);
%     error_lap = max(abs(lapH_true - lapH((0*m*n+1001):(0*m*n+1050))))/max(abs(lapH_true))
%     n_true = reshape(n_true, [50, 3]);
%     fb_true = repmat(lapH_true + 2*H_true.*(H_true.^2 - K_true), [1, 3]).*n_true
%     error_fb = max(vecnorm(fb_true - fb((0*m*n+1001):(0*m*n+1050),:),2,2))/max(vecnorm(fb_true,2,2)), 
%     
    H = -H;
    lapH = S.getLaplacian(H);
    
    max_H = max(H)
    max_H2 = max(H.^2)
    max_H2K = max(H.^2-K)
    max_K= max(K)
    maxhk = max(H.*(H.^2-K))
    maxlapH = max(lapH)
    maxlapHplus = max(lapH + 2*H.*(H.^2-K))
    max_norm_fb = max(vecnorm(-fb,2,2))
    %all_r = [];
    for ii=1:6
        
        S.patches(ii).q_sl = fb((ii-1)*m*n+1:ii*m*n, :);
        %all_r = [all_r; S.patches(ii).r];
        
    end
%     all_r = all_r';
%     all_r = all_r(:);
%     dlmwrite('ellipsoid2_matlabtrg32.txt',all_r,'newline','pc', 'precision', 12);
%     figure;
%     for i=1:6
%         p = patches(i);
%         scatter3(p.r(:,1),p.r(:,2),p.r(:,3));
%         hold on;
%     
%     end
    SLerror = zeros(size(all_trg, 1), 1);
    true_vals = zeros(size(all_trg, 1), 1);
    true_norm = max(vecnorm(all_pot, 2, 2));
    for ii=1:size(all_trg, 1)
        trg = all_trg(ii, :);
        true_val = all_pot(ii, :);
        norm_trg = sqrt(norm(trg)) ;
        val = SLSmooth(trg,  S);
        pot = all_pot(ii, :);
        SLerror(ii) = norm(val + true_val);
        true_vals(ii) = norm(true_val);
        
    end
   
    error = max(SLerror)/true_norm
end

