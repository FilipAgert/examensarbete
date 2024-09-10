function coord = getCoord(idx, cols)
    row = floor(idx/cols-1e-10)+ 1;
    col = mod(idx, cols);
    if col == 0
        col = cols;
    end
    coord = [row, col];
end