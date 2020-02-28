function f_fft = ves2fft(f, nu, nv)
%Function to convert vesicle code vector values to fft format

        f_mat = reshape(f, [nu nv]); %u varies column wise, v rowwise
        f_mat = [zeros(1, nv);f_mat];
        f_mat = [zeros(nu+1, 1) f_mat];
        f_fft = 2*f_mat'; %2 for scaling due to transformation [0, pi]->[-pi, pi]
    



end