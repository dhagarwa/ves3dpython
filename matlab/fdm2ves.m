function f_ves = fdm2ves(f, nu, nv)
%Function to convert fdm vector values to ves format

        f_mat = reshape(f, [nv nu]); %u varies column wise, v rowwise
        f_ves = f_mat'; %no scaling here, we do it directly while extracting f_du, f_dv 
        f_ves = f_ves(:);
    

end