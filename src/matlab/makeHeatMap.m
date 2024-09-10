function makeHeatMap(weights, cols)
    A = zeros(length(weights)/cols, cols);
    for i = 1 : length(weights)
        coord = getCoord(i, cols);
        A(coord(1), coord(2)) = weights(i);
    end
    colormap('hot');
    imagesc(A);
    colorbar;
end
