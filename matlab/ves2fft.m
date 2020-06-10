function f_fft = ves2fft(f, nu, nv)
%Function to convert vesicle code vector values to fft format

        f_mat = reshape(f, [nu nv]); %u varies column wise, v rowwise
        
        f_mat = [zeros(1, nv);f_mat];
        f_mat = [zeros(nu+1, 1) f_mat];
        f_fft = f_mat'; %no scaling here, we do it directly while extracting f_du, f_dv 
    



end