function F_filtered = gaussian_filter(F)
% f is 2D matrix of spectrum of a 2d function 
% apply 2D gaussian high pass filter on it 
    
    m = size(F, 1);
    n = size(F, 2);
    H = zeros(m, n);
    D0 = (m)/2+2; %cutoff frequency
    k0 = m/2;
    l0 = n/2;
    a = 1; %attenuating factor 
    for ii = 1:m        
        for jj = 1:n            
            %H(ii, jj) = 1 - exp(-a*((ii-k0)^2+(jj-l0)^2)/D0^2);
            if (ii-k0)^2 + (jj-l0)^2 > D0^2 
                H(ii, jj) = 1;
            end
        end
        
    end
    
    F_filtered = H.*F;


end