close all;
clear all;
1;
file0 = '../../data/RUN-4';
file1 = '../../data/RUN-15';
fileActual = '../../data/RUN-10';
dat0 = readmatrix(file0, 'Delimiter', ';');
dat1 = readmatrix(file1, 'Delimiter', ';');
datActual = readmatrix(fileActual, 'Delimiter', ';');

Ncoords = size(dat0,1)*size(dat0,2);
x = zeros(2,Ncoords);
pdf0 = zeros(1,Ncoords);
pdf1 = zeros(1,Ncoords);
pdf2 = zeros(1,Ncoords);
idx = 1;
for i = 1:size(dat0,1)
    for j = 1:size(dat1,2)
        x(1,idx) = i;
        x(2,idx) = j;%allocate coordinates.
        pdf0(idx) = dat0(i,j);
        pdf1(idx) = dat1(i,j);
        pdf2(idx) = datActual(i,j);
        idx = idx + 1;
    end
end

s = eye(2);
E0 = 4/21;
E1 = 15/21;
E = 12/21;
pdf0 = pdf0/sum(pdf0);%assure normalized
pdf1 = pdf1/sum(pdf1);%assure normalized
pdf2 = pdf2/sum(pdf2);
subplot(1,5,1)
plotheatmap(x,pdf0);
subplot(1,5,2)
plotheatmap(x,pdf1);
subplot(1,5,3)
plotheatmap(x,pdf2);
subplot(1,5,4)

pdFast = guessPdfFast(x,E,E0,E1,pdf0,pdf1);
plotheatmap(x,pdFast);
subplot(1,5,5)
pd = guessPdf(x,E,E0,E1,pdf0,pdf1);
plotheatmap(x,pd);
%scatter3(x(1,:), x(2,:), pdf0);
%hold on;
%scatter3(x(1,:), x(2,:), pdf1);
%scatter3(x(1,:), x(2,:), pd);
%%
close all;
clear all;
1;
file0 = '../../data/RUN-4';
file1 = '../../data/RUN-15';
fileActual = '../../data/RUN-10';
dat0 = readmatrix(file0, 'Delimiter', ';');
dat1 = readmatrix(file1, 'Delimiter', ';');
datActual = readmatrix(fileActual, 'Delimiter', ';');
pdf0 = dat0(51,:);
pdf1 = dat1(51,:);
pdf2 = datActual(51,:);

pdf0=pdf0/sum(pdf0);
pdf1=pdf1/sum(pdf1);
pdf2=pdf2/sum(pdf2);
x = 1:length(pdf0);
%pdf0 = pdf('Normal', x,20,5);
%pdf1 = pdf('Normal', x, 60,5);
%pdf2 = pdf('Normal', x, 40, 5);
pdf0=pdf0/sum(pdf0);
pdf1=pdf1/sum(pdf1);
pdf2=pdf2/sum(pdf2);



Ncoords = size(dat0,1)*size(dat0,2);

E0 = 2
E1 = 6;
E = 4;
plot(x,pdf0, 'b--');
hold on;
plot(x,pdf1, 'b--');
plot(x,pdf2, 'r--');
pdFast = guessPdfFast(x,E,E0,E1,pdf0,pdf1);
pd = guessPdf(x,E,E0,E1,pdf0,pdf1);
plot(x,pd, 'g')
plot(x,pdFast, 'k')

function plotheatmap(X, values)
    % Assuming you have the vectors x, y, and values
    x = X(1,:); % Example x-coordinates
    y = X(2,:); % Example y-coordinates
    numbins = sqrt(length(values));
    % Define the grid resolution (e.g., 100x100 grid)
    x_edges = linspace(min(x), max(x), 0.5*round(numbins)); % x-axis grid
    y_edges = linspace(min(y), max(y), 0.5*round(numbins)); % y-axis grid

    % Bin the data into the grid (using histcounts2 for counting and averaging)
    [~, ~, x_bin] = histcounts(x, x_edges);
    [~, ~, y_bin] = histcounts(y, y_edges);

    % Initialize the grid for heatmap values
    heatmap_grid = nan(length(y_edges)-1, length(x_edges)-1);
    counts = zeros(size(heatmap_grid));

    % Loop through each point and add the value to the corresponding grid cell
    for i = 1:length(values)
        if x_bin(i) > 0 && y_bin(i) > 0  % Ensure the bin index is valid
            heatmap_grid(y_bin(i), x_bin(i)) = nansum([heatmap_grid(y_bin(i), x_bin(i)), values(i)]);
            counts(y_bin(i), x_bin(i)) = counts(y_bin(i), x_bin(i)) + 1;
        end
    end

    % Average the values for each grid cell
    heatmap_grid = heatmap_grid ./ counts;

    % Plot the heatmap
    imagesc(x_edges, y_edges, heatmap_grid);
    axis xy; % Ensure the axis is correctly oriented
    colorbar; % Display color bar to show value scale
    title('Heatmap of values over 2D space');
    xlabel('X-axis');
    ylabel('Y-axis')

end
function GUESS = guessPdfFast(x,E,E0,E1,pdf0,pdf1) %can use when E close to E0, E1.
    a = alpha(E, [E0, E1]);
    GUESS = (1-a)*pdf0 + a*pdf1;
    GUESS = GUESS/sum(GUESS);
end
function GUESS = guessPdf(x, E, E0, E1, pdf0, pdf1)
%X Is N-dimensional coordinate of size NxM where M is number of grid points
%E0, E1 are the energies of the probability densities
%E is energy of probability density you want to evaluate
%pd0, pd1 are 1xM probability densities associated to each point in X.
    
    a = alpha(E, [E0, E1]);
    mu0 = getMean(x,pdf0);
    K0 = getCovariance(x,pdf0);
    mu1 = getMean(x, pdf1);
    K1 = getCovariance(x,pdf1);
    mu = (1-a)*mu0 + a*mu1;
    K = (1-a)*K0 + a*K1;
    
    deltaMu = abs(mu0-mu1);
    muMult = deltaMu*deltaMu';
    Kval = det(K);
    disp("muval: " + max(muMult));
    disp("K val: " + Kval);
    
    
    
    K0Kinv = K0/K;
    K1Kinv = K1/K;
    x0 = sqrtm(K0Kinv)*(x-mu) + mu0; %K0/K = K0 * inv(K)
    x1 = sqrtm(K1Kinv)*(x-mu) + mu1;
    pd = 0;
    gamma0 = genGammaMatrix(x,x0,1);
    gamma1 = genGammaMatrix(x,x1,1);
    for i = 1:size(x,2)
        pd = pd + (1-a)*gamma0(i,:)*pdf0(i) + a*gamma1(i,:)*pdf1(i);
    end
    pd = pd/sum(pd);
    GUESS = pd;

end


function a = alpha(inputCoordinate, coords)
    a = (inputCoordinate-coords(1))/(coords(2)-coords(1));
end 

function mu = getMean(x, val)
    %Let val be 1-dimensional of length M
    %Let x be NxM dimensional
    %Where N is number of dimensions.
    mu = 0;
    for i = 1:size(x, 2)
        coord = x(:,i);
        mu = mu + coord*val(i);
    end
end

function c = getCovariance(x, val)
    %Let val be 1-dimensional of length M
    %Let x be NxM dimensional
    %Where N is number of dimensions.
    c = 0;
    mu = getMean(x,val);
    
    for i = 1:size(x,2)
        coord = x(:,i);
        
        u =  coord-mu;
        c = c + u*u'*val(i);
    end
end

%TODO: generate nearest neighbour for point.
function gammaMat = genGammaMatrix(x, x1, nSamples) %O(n^2)
    %Define the hypercube for X
    %Sample length(x) * 10 times
    %For each sample, check nearest neighbour in x and in x1
    %add to matrix.
    numDims = size(x,1);
    
    
    rangeX = [min(x,[],2),max(x,[],2)];
    rangeX1 = [min(x1,[],2),max(x1,[],2)];
    gammaMat = zeros(length(x), length(x1));
    kSpace = nSamples*size(x,2);
    for k = 1:(kSpace)
        coord = zeros(numDims,1);
        for dim = 1:numDims
            coord(dim) = randbetween(rangeX1(dim,1), rangeX1(dim,2));
        end
        
        nearestInX = nearestNeighbour(coord, x);
        nearestInX1 = nearestNeighbour(coord, x1);
        gammaMat(nearestInX, nearestInX1) = gammaMat(nearestInX,nearestInX1) + 1;
        if(mod(k,round(kSpace/10)) == 0)
            disp(k + " out of " + kSpace)
        end
    end
    s = sum(gammaMat,1);
    zeroSum = find(s==0);
    for idx = zeroSum
        coord = x1(:,idx);
        nearestInX = nearestNeighbour(coord, x);
        gammaMat(nearestInX, idx) = gammaMat(nearestInX,nearestInX1) + 1;
    end
    s = sum(gammaMat,1);
    gammaMat = gammaMat./s;
    
end




function r = randbetween(a,b)
    r = (b-a).*rand + a;

end
function d = dist(a,b, fullDist)
    %d = norm((a-b)./fullDist);
    d = norm(a-b);
end

function [idx] = nearestNeighbour(p, ps)
    dis = inf;
    fullDist = max(ps, [], 2)-min(ps,[],2);
    for i = 1:length(ps)
        neighbour = ps(:,i);
        currentDis = dist(p, neighbour, fullDist);
        if(currentDis < dis)
            dis = currentDis;
            idx = i;
            
        end
    end
end
