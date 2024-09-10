function A = linkCells(A, fromCoord, toCoord, rows, cols)
%A is probability matrix.
    fromIndex = getIndex(fromCoord, cols);
    toIndex = getIndex(toCoord, cols);
    
    neighbours = find(A(:,fromIndex));
    newProb = 1/(length(neighbours) + 1);
    A([neighbours; toIndex], fromIndex) = newProb;
    
end

