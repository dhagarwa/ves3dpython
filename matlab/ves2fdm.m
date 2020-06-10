function f_fdm = ves2fdm(f, nu, nv)
%Function to convert vesicle code vector values to fdm format

        f_mat = reshape(f, [nu nv]); %u varies column wise, v rowwise
        f_fdm = f_mat'; %no scaling here, we do it directly while extracting f_du, f_dv 
        f_fdm = f_fdm(:);
    

end