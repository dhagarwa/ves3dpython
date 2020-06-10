function [f_du, f_dv] = derivFDM(f, NU, NV, h_u, h_v, ood, ooa)

    [Du, Dv] = getNonCompactFDmatrix2D(NU, NV, h_u, h_v, ood,ooa);
    % apply the diff. matrix
    f_du = Du*f;
    f_dv = Dv*f;


end