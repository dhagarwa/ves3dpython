function f_interpolated = upsample(f, patch)
% upsample a periodic function f values on 2D standard patch

    %uf = patch.uf; %upsampling factor 
    
    %m = patch.Nu;
    %n = patch.Nv;
    
    uf = 2; m = 15; n = 15;
    
    f_mat = reshape(f, [m, n]);
    f_mat = [zeros(1, n);f_mat];
    f_mat = [zeros(m+1, 1) f_mat];
    stem3(f_mat);
    figure;
    stem3(abs(fft2(f_mat)));
    f_mat_colpad = [];
    %zero pad the f_mat for upsampling
    for jj=1:n+1
       f_mat_colpad = [f_mat_colpad  f_mat(:, jj) zeros(m+1, uf-1)];
                
    end
    %f_mat_colpad = [f_mat_colpad zeros(m, 1)];
    
    f_mat_pad = [];
    for ii=1:m+1
        f_mat_pad = [f_mat_pad;  f_mat_colpad(ii, :); zeros(uf-1, uf*(n+1))];        
    end
    
    %f_mat_pad = [f_mat_pad; zeros(1, uf*n+1)];
    
    size(f_mat_pad) 
    figure;
    stem3(f_mat_pad);
    %get fourier transform of padded f
    F_mat_pad = fft2(f_mat_pad);
    
    %not centralizing the spectrum F_mat_pad so using high pass filter. if centralizing use low
    %pass filter. 
    figure;
    stem3(abs(F_mat_pad));
    %f_mat_interpolated_nofilter = real(ifft2(F_mat_pad));
    %figure;
    %stem3(f_mat_interpolated_nofilter);
    
    F_mat_pad = gaussian_filter(F_mat_pad); 
    figure;
    stem3(abs(F_mat_pad));
    %once filtered take ifft
    
    f_mat_interpolated = real(ifft2(F_mat_pad));
    
    f_interpolated = f_mat_interpolated(:);
    
    figure;
    stem3(f_mat_interpolated);



end