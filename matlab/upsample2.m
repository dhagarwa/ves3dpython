function f_interpolated = upsample2(f, patch)
% upsample a periodic function f values on 2D standard patch

    %uf = patch.uf; %upsampling factor 
    
    %m = patch.Nu;
    %n = patch.Nv;
    
    uf = patch.uf; m = patch.Nu; n = patch.Nv;
    
    f_mat = reshape(f, [m, n]);
    f_mat = [zeros(1, n);f_mat];
    f_mat = [zeros(m+1, 1) f_mat];
    %figure;
    %stem3(f_mat);
    
    
    %take fft2 of this f_mat
    
    F_mat = fft2(f_mat); %2D dft of f_mat
    %figure;
    %stem3(abs(F_mat));
    %stuff fft with zeros in high frequency region
    
    z = zeros((m+1)/2, (n+1)/2);
    
    F_mat_padded = zeros(uf*(m+1), uf*(n+1));
    %F_mat_padded = [F_mat(1:(m+1)/2, 1:(n+1)/2) 0.5*F_mat(1:(m+1)/2, (n+1)/2+1) z(:, 2:end) z 0.5*F_mat(1:(m+1)/2, (n+1)/2+1) F_mat(1:(m+1)/2, (n+5)/2:n+1); 0.5*F_mat((m+3)/2, 1:(n+3)/2) z(1, 2:end) z(1, :) 0.5*F_mat((m+3)/2, ); z(2:end, :) z(2:end, :) z(2:end, :) z(2:end, :); z z z z; F_mat((m+3)/2:m+1, 1:(n+1)/2) z z F_mat((m+3)/2:m+1, (n+3)/2:n+1)];
    F_mat_padded = [F_mat(1:(m+1)/2, 1:(n+1)/2) z z F_mat(1:(m+1)/2, (n+3)/2:n+1); z z z z; z z z z; F_mat((m+3)/2:m+1, 1:(n+1)/2) z z F_mat((m+3)/2:m+1, (n+3)/2:n+1)];
    
    
    jj = (n+1)/2+1;
    for ii=1:uf*(m+1)
        if ii == (m+1)/2+1            
            F_mat_padded(ii, jj) = F_mat((m+1)/2+1,jj)/4;
                        
        end
       
        if ii<(m+1)/2+1
            F_mat_padded(ii, jj) = F_mat(ii,jj)/2; 
            
        end
         
        if ii==3*(m+1)/2+1
            F_mat_padded(ii, jj) = F_mat((m+1)/2+1,jj)/4;
        end
        
        if ii>3*(m+1)/2+1
           F_mat_padded(ii, jj) = F_mat(ii-3*(m+1)/2+(m+1)/2, jj)/2;
            
        end
        
    end

    jj = 3*(n+1)/2+1;
    for ii=1:uf*(m+1)
        if ii == (m+1)/2+1            
            F_mat_padded(ii, jj) = F_mat((m+1)/2+1,(n+1)/2+1)/4;
                        
        end
       
        if ii<(m+1)/2+1
            F_mat_padded(ii, jj) = F_mat(ii,(n+1)/2+1)/2; 
            
        end
         
        if ii==3*(m+1)/2+1
            F_mat_padded(ii, jj) = F_mat((m+1)/2+1,(n+1)/2+1)/4;
        end
        
        if ii>3*(m+1)/2+1
           F_mat_padded(ii, jj) = F_mat(ii-3*(m+1)/2+(m+1)/2, (n+1)/2+1)/2;
            
        end
        
    end
    
    ii = (m+1)/2+1;
    for jj=1:uf*(n+1)

       
        if jj<(n+1)/2+1
            F_mat_padded(ii, jj) = F_mat(ii,jj)/2; 
            
        end
         
        
        if jj>3*(n+1)/2+1
           F_mat_padded(ii, jj) = F_mat(ii, jj-3*(n+1)/2+(n+1)/2)/2;
            
        end
        
    end

    ii = 3*(m+1)/2+1;
    for jj=1:uf*(n+1)
        
        if jj<(n+1)/2+1
            F_mat_padded(ii, jj) = F_mat((m+1)/2+1,jj)/2; 
            
        end
         
        
        if jj>3*(n+1)/2+1
           F_mat_padded(ii, jj) = F_mat((m+1)/2+1, jj-3*(n+1)/2+(n+1)/2)/2;
            
        end
        
    end
    
    
    %F_mat_padded(1:8,1:8) = F_mat(1:8,1:8);
    %F_mat_padded(3*(n+1)/2+2:end,1:8) = F_mat((n+5)/2:end,1:8);
    
    %F_mat_padded;
    %figure; 
    %stem3(abs(F_mat_padded));
    f_interpolated = real(ifft2(F_mat_padded))*uf*uf;
    f_interpolated = f_interpolated(2:end, 2:end);
    %figure;
    %stem3(f_interpolated);
    f_interpolated = f_interpolated(:);
    
end