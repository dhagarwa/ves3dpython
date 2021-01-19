function [] = cylinder_theta_phi()
    %Function to generate the cylinder surface
    %
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'theta_out16.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    u = cell2mat(C);
    
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'phi_out16.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    v = cell2mat(C);
    
    
    %u = u(:);
    %v = v(:);
    r = zeros(size(u,1),3);
    for ii=1:size(u, 1)
       u0 = u(ii); v0 = v(ii);
       if ii <= 128
           r(ii,1) = -8;
           r(ii,2) = (floor(ii/32)+1)*5*sin(v0);
           r(ii,3) = (floor(ii/32)+1)*5*sin(v0);
           
       elseif ii>13*128
           r(ii,1) = 12;
           r(ii,2) = (floor(ii/32)+1)*5*sin(v0);
           r(ii,3) = (floor(ii/32)+1)*5*sin(v0);
              
        
    end
    scatter3(r(:,1),r(:,2),r(:,3));
    r = r(:); %C++ src corrdinate form
    
    %write r to text file
    dlmwrite('cylinder_src16.txt',r,'newline','pc');
end
