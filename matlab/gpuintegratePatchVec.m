function val = gpuintegratePatchVec(p, f)

    val_x = gpuintegratePatch(p, f(:, 1));
    val_y = gpuintegratePatch(p, f(:, 2));
    val_z = gpuintegratePatch(p, f(:, 3));
    val = [val_x, val_y, val_z];


end
