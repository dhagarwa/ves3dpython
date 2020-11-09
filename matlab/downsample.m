function f_down = downsample(f_in, patch)
% downsample a periodic function f values on 2D standard patch
% f is given in upsampled fft matrix form 
%transform it to fft matrix form but downsampled
    %uf = patch.uf; %upsampling factor 
    
    %m = patch.Nu;
    %n = patch.Nv;
    
    uf = patch.uf; m = patch.Nu; n = patch.Nv;
    
    
    mid = (m+1)/2 +1;
    next_mid = 3*(mid-1) +1;
    f = fft2(f_in);
    f_down = [f(1:mid,1:mid) f(1:mid,next_mid+1:end); f(next_mid+1:end,1:mid) f(next_mid+1:end, next_mid+1:end)];
    
    ii = mid;
    for jj=1:n+1
       
       f_down(ii, jj) = f_down(ii, jj)*2;
       if jj == mid
           f_down(ii, jj) = f_down(ii, jj)*2;
       end
                 
    end
    
    jj = mid;
    
    for ii=1:m+1
       
     
       if ii ~= mid
           f_down(ii, jj) = f_down(ii, jj)*2;
       end
                 
    end    
    
    f_down = real(ifft2(f_down))/(uf*uf);
    
end