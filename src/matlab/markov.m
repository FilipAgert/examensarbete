clear;
close all;
rows = 11;
cols = 11;
A = generateMatrix(rows, cols);
startWaveFunction = zeros(rows*cols,1);
startcell = getMiddleCell(rows, cols);
startIndex = getIndex(startcell, cols);
startWaveFunction(startIndex) = 1;
    
weights = startWaveFunction;
 
A = linkCells(A, [5, 1], getMiddleCell(rows, cols), rows, cols); 
A = linkCells(A, [6, 1], getMiddleCell(rows, cols), rows, cols); 
A = linkCells(A, [7, 1], getMiddleCell(rows, cols), rows, cols); 

A = linkCells(A, [6, cols], getMiddleCell(rows, cols), rows, cols); 

mat = A;
for t = 1 : 500
    mat = A * mat;
end

for t = 1 : 5001
    makeHeatMap(weights, cols);
    if t == 1
        pause (5)
    end
    weights = A * weights;
    pause(0.01);
end

function A = generateMatrix(rows, cols)
    for r = 1 : rows
       for c = 1 : cols
           coord = [r, c];
           idx = getIndex(coord, cols);
           adjecent = getNeighbours(coord, rows, cols, false);
           adjIdxs = [];
           for adjIdx = 1 : size(adjecent,1)
               adj = adjecent(adjIdx,:);
               adjecentIndex = getIndex(adj,cols);
               adjIdxs = [adjIdxs, adjecentIndex];
           end
           prob = getEqualProbabilities(length(adjIdxs));
           for n = 1 : length(adjIdxs)
               A(adjIdxs(n),idx) = prob(n);
           end
       end
    end
end





function coords = getNeighbours(coord, rows, cols, includeCorners)
    rCords = [coord(1) - 1, coord(1), coord(1) + 1];
    cCords = [coord(2) - 1, coord(2), coord(2) + 1];
    [p, q] = meshgrid(rCords, cCords);
    coords = [p(:), q(:)];
    if(~includeCorners)
        coords(1,:) = [];
        coords(2,:) = [];
        coords(end,:) = [];
        coords(end - 1, :) = [];
    end
    [belowLowRows, ~] = find(coords < 1);
    [row2, col2] = find(coords > rows);
    relevantIdx = find(col2 == 1);
    row2 = row2(relevantIdx);
    [row3, col3] = find(coords > cols);
    relevantIdx = find(col3 == 2);
    row3 = row3(relevantIdx);
    removedRows = [belowLowRows; row2; row3];
    coords(removedRows,:) = [];
end

function probs = getRandomProbabilities(number)
    probs = rand(1,number);
    probs = probs/sum(probs);
end

function probs = getEqualProbabilities(number)
    probs = ones(1,number);
    probs = probs/number;
end