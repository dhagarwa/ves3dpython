function f_ves = fft2ves(f, nu, nv)
%Function to convert f from fft matrix form to vesicle code column
%form
    
        f_ves = f';
        f_ves = f_ves(2:nu, 2:nv);
        f_ves = f_ves(:);
        


end