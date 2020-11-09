function [] = flower_surface()
    %Function to generate the flower surface
    % r = [psin(u)cos(v);psin(u)sin(v);pcos(u)]
    % p = 1 + exp(-3 Real(Y(3,2,u,v))), u \in [0,pi], v \in [0,2*pi]
    m = 32;
    n = 2*32;
    [u,v]=meshgrid(pi*(0:m)/(m+1),  (2*pi)*(0:n-1)/(n)); %u and v both are m x n matrices

    u = u(:);
    v = v(:);
    
    r = flower_function(u,v);
    scatter3(r(:,1),r(:,2),r(:,3));
    r = r(:); %C++ src corrdinate form
    
    %write r to text file
    dlmwrite('flower_src.txt',r,'newline','pc');
end

