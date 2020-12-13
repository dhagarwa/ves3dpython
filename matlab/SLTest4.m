function [] = SLTest4()
      clc

    m = 31;
    n = 31;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'flower_trg.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    A = cell2mat(C);
    m0 = size(A, 1);
    all_trg = reshape(A, [3, m0/3]);
    all_trg = all_trg';

    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'matlab_flower_SL_potential_trg.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    A = cell2mat(C);
    m0 = size(A, 1);
    all_pot = reshape(A, [3, m0/3]);
    all_pot = all_pot';
    patches = [];
    trg = [0.1690   -0.0903   -2.0063]; val = 0;
    for i=1:6
       patch =  flowerPatch(m, n, i, R);
       
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
    S = Surface(patches, [0, 0, 0]);
    DLvalue = DLSmooth(trg,1,S)
%     figure;
%     for i=1:6
%         p = patches(i);
%         scatter3(p.r(:,1),p.r(:,2),p.r(:,3));
%         hold on;
%     
%     end
    SLerror = zeros(size(all_trg, 1), 1);
    true_vals = zeros(size(all_trg, 1), 1);
    for ii=1:size(all_trg, 1)
        trg = all_trg(ii, :);
        true_val = all_pot(ii, :);
        norm_trg = sqrt(norm(trg)) ;
        val = SLSmooth(trg,  S);
        pot = all_pot(ii, :);
        SLerror(ii) = norm(val - true_val);
        true_vals(ii) = norm(true_val);
        
    end
   
    error = max(SLerror)/max(true_vals)
end

