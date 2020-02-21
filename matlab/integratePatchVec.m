function val = integratePatchVec(p, f)

    val_x = integratePatch(p, f(:, 1));
    val_y = integratePatch(p, f(:, 2));
    val_z = integratePatch(p, f(:, 3));
    val = [val_x, val_y, val_z];


end