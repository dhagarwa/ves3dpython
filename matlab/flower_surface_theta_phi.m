function [] = flower_surface_theta_phi()
    %Function to generate the flower surface
    % r = [psin(u)cos(v);psin(u)sin(v);pcos(u)]
    % p = 1 + exp(-3 Real(Y(3,2,u,v))), u \in [0,pi], v \in [0,2*pi]
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'theta_out48.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    u = cell2mat(C);
    
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'phi_out48.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    v = cell2mat(C);
    
    
    %u = u(:);
    %v = v(:);
    
    r = flower_function(u,v);
    scatter3(r(:,1),r(:,2),r(:,3));
    r = r(:); %C++ src corrdinate form
    
    %write r to text file
    dlmwrite('flower_src48.txt',r,'newline','pc');
end

