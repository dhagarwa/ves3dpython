function f = stripboundary(f_, p)
%strip boundary values and return column vector of values of f at interior
%nodes only in patch p 

    m = p.Nu; n = p.Nv;
    
    f_mat = reshape(f_,[m+2,n+2]);
    f_mat = f_mat(2:m+1,2:n+1);
    f = f_mat(:);
    


end