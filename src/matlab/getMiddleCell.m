function coord = getMiddleCell(row, col) %Gets center cell of matrix Row x col
    coord = [floor((row - 1)/2 + 1), floor((col - 1)/2 + 1)];
end
