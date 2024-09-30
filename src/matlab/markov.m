clear;
close all;
rows = 40;
cols = 40;
n = rows*cols;
tol = 1/(n*1e3);
steps = 0;
A = generateMatrix(rows, cols);
startWaveFunction = zeros(rows*cols,1);
startcell = getMiddleCell(rows, cols);
startIndex = getIndex(startcell, cols);
startWaveFunction(startIndex) = 1;
    
weights = startWaveFunction;
 
A = linkCells(A, [round(rows/2-1), 1], getMiddleCell(rows, cols), rows, cols); 
A = linkCells(A, [round(rows/2), 1], getMiddleCell(rows, cols), rows, cols); 
A = linkCells(A, [round(rows/2+1), 1], getMiddleCell(rows, cols), rows, cols); 

A = linkCells(A, [round(rows/2-1), cols], getMiddleCell(rows, cols), rows, cols); 
A = linkCells(A, [round(rows/2), cols], getMiddleCell(rows, cols), rows, cols); 
A = linkCells(A, [round(rows/2+1), cols], getMiddleCell(rows, cols), rows, cols); 
converged = false;
T1 = datetime;
while ~(converged)
    prev = weights;
    weights = A * weights;
    steps = steps + 1;
    if(mod(steps,10) == 0)
        if(norm(prev-weights)<tol)
            converged = true;
        end
    end
    if(mod(steps, 100) == 0)
        disp(steps)
    end
    
end
T2 = datetime;
B = (A - eye(n,n));
sol = null(B, tol*1e2);
sol = sol./sum(sol);
T3 = datetime;
zer = zeros(n,1);
restart = 100;
zer(1) = 1;
tolGMRES = 0.015;

a21 = B(n,1:(end-1));
a12 = B(1:(end-1),n);
C = B(1:(end-1), 1:(end-1));
D  = sparse(C);
[L, U] = ilu(D, struct('type','ilutp','droptol', 1e-6));
tempSol = gmres(D, -a12, [], tolGMRES, 500, L,U, ones(n-1,1)/n);
solPart = [tempSol; 1];
solPart = solPart./sum(solPart);
T4 = datetime;
bicgSol = bicg(D,-a12, tolGMRES, 100,L, U, ones(n-1,1)/n);
bicgFullSol = [bicgSol ; 1];
bicgFullSol = bicgFullSol./sum(bicgFullSol);
T5 = datetime;



subplot(3,2,1)
title('linsolve')
makeHeatMap(sol, cols)
subplot(3,2,2)
title('iterative')
makeHeatMap(weights, cols)
subplot(3,2,3)

makeHeatMap(solPart, cols)
subplot(3,2,4)
makeHeatMap(bicgFullSol, cols)
IterativeTime = seconds(T2-T1)
Linsystemtime = seconds(T3-T2)
gmrestime = seconds(T4-T3)
bicgtime = seconds(T5-T4)

function A = generateMatrix(rows, cols)
    n = rows*cols;
    nnz = 5*n;
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