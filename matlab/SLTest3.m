function [] = SLTest3()
      clc

    m = 63;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'matlab_test_SL_targets.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    A = cell2mat(C);
    m0 = size(A, 1);
    all_trg = reshape(A, [m0/3, 3]);

    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'matlab_test_SL_potential.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    A = cell2mat(C);
    m0 = size(A, 1);
    all_pot = reshape(A, [m0/3, 3]);
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S = Surface(patches, [0, 0, 0]);
 
    SLerror = zeros(size(all_trg, 1), 1);
    for ii=1:60:size(all_trg, 1)
        trg = all_trg(ii, :);
        true_val = all_pot(ii, :);
        norm_trg = sqrt(norm(trg)) ;
        val = SLSmooth(trg,  S);
        pot = all_pot(ii, :);
        SLerror(ii) = norm(val - true_val)/norm(true_val);
        
    end
   
    max(SLerror)
end

