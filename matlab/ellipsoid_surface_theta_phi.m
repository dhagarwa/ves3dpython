function [] = ellipsoid_surface_theta_phi()
    %Function to generate the ellipsoid surface

    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'theta_out64.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    u = cell2mat(C);
    
    rootdir = '/Users/dhwanit/Google Drive/Biros Research/ves3dpython/matlab';
    relativefolder = 'testfiles';
    fid = fopen(fullfile(rootdir, relativefolder, 'phi_out64.txt'), 'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    %celldisp(C);
    v = cell2mat(C);
    
    
    %u = u(:);
    %v = v(:);
    a = 1; b = 0.9; c = 1;
    r = ellipsoid_function(a, b, c, u,v);
    scatter3(r(:,1),r(:,2),r(:,3));
    
    %r = r'; r = r(:);
    %dlmwrite('ellipsoid2_trg64.txt',r,'newline','pc', 'precision', 16);
    
    r = r(:); %C++ src corrdinate form
    
    %write r to text file
    dlmwrite('ellipsoid3_src64.txt',r,'newline','pc', 'precision', 16);
end

function r = ellipsoid_function(a, b, c, u, v)
    z = c*cos(u);
    x = a*sin(u).*cos(v);
    y = b*sin(u).*sin(v);
    r = [x y z];

end

