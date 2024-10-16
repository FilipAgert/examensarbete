clear

figure(11)
MIN_II = 1;
MAX_II = 53;
MIN_JJ = 1;
MAX_JJ = 15;
MIN_KK = 1;
MAX_KK = 15;
MIN_LL = 1;
MAX_LL = 15;
MIN_MM = -40;
MAX_MM = 40;

%MIN_MAX_DIM = [MIN_II, MAX_II; MIN_JJ, MAX_JJ; MIN_KK, MAX_KK; MIN_LL, MAX_LL; 0, MAX_MM];
%dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ, 1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM];
MIN_MAX_DIM = [MIN_II, MAX_II; MIN_JJ, MAX_JJ; MIN_KK, MAX_KK; MIN_LL, MAX_LL; MIN_MM, MAX_MM];
dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ, 1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM - MIN_MM];

startCoord = [17,7,3,13,30];
%pd = str2double(pdstr)
%This is very wrong.
AZ = 102;
AA = 256;

A = zeros(dimSize);
files = ["../../data/PDMM-20","../../data/PDMM-25", "../../data/PDMM-30", "../../data/PDMM-40"];
fusionChance = [39.5, 80.1,58.9,46.4];
startingCoords = [19, 6, 3, 13, 30; 
                  18, 7, 3, 14, 29; 
                  17, 7, 3, 13, 30;
                  16, 7, 3, 13, 30];
E = [20,25,30,40];
for f = 1:4
    subplot(2,2,f)
    pd = readmatrix(files(f), "Delimiter", ";");
    
    
    
    for i = 1:length(pd)
        if(mod(i,1000000) == 0)
            disp(100*i/length(pd) + " percent complete")
        end
        c = coordFromLinearIdxShifted(i, dimSize, MIN_MAX_DIM);
        A(c(1),c(2),c(3),c(4),c(5)) = pd(i);
    end
    A(find(A<0)) = 0;
    
    maxval = 0;
    s = 0;
    maxj = 0;
    maxk = 0;
    maxl = 0;
    summedA = zeros([dimSize(1), dimSize(5)]);
    TEMP = permute(A,[1,5,2,3,4]); %Move dim5 into dim 2. For plotting purposes I, M, J, K, L
    for j = MIN_JJ:MAX_JJ
        for k = MIN_KK:MAX_KK
            for l = MIN_LL:MAX_LL
                summedA(:,:) = summedA(:,:) + TEMP(:,:,j,k,l);
            end
        end
    end
    
    xrang = (MIN_MAX_DIM(5,1):MIN_MAX_DIM(5,2))*0.02;
    yrang = MIN_MAX_DIM(1,1):MIN_MAX_DIM(1,2);
    for i = 1:length(yrang)
        val = get_q2_val(AZ,AA,yrang(i));
        yrang(i) = val;
    end
    %Now we have to create a kernel
    kernel = zeros([dimSize(1), dimSize(5)]);
    [X,Y] = meshgrid(xrang, yrang);
    bw = 0.003
    bwx = bw*(xrang(end)-xrang(1));
    bwy = bw*(yrang(end)-yrang(1));
    summedA(find(summedA<1e-12))=1e-12;
    for i = 1:MIN_MAX_DIM(1,2)
        for m = MIN_MAX_DIM(5,1):MIN_MAX_DIM(5,2)
            ms = m - MIN_MM + 1;
            p = summedA(i,ms);
            y = yrang(i);
            x = xrang(ms);

            kernel = kernel + p*exp(-((X-x).^2)/(2*bwx^2) -( (Y-y).^2)/(2*bwy^2));

        end
    end



    
    plotMatrix = log10(summedA);
    %plotMatrix = kernel;
    climits = [min(plotMatrix,[],'all'), max(plotMatrix,[],'all')];
    
    
    x_coord = startingCoords(f,5)*0.02;
    y_coord = get_q2_val(AZ,AA,startingCoords(f,1));


    plot(x_coord, y_coord, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'm');
    hold on
    pcolor(X,Y,plotMatrix);
    plot(x_coord, y_coord, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'm');
    shading interp;
    
    %imagesc('XData', xrang, 'YData', yrang,'CData', plotMatrix);
    %contourf(xrang, yrang, plotMatrix, 20)
    
    
    
    
    
    %legend('Probability density', 'Starting coordinate')
    legend('Starting coordinate')
    clim(climits);
    xlim([MIN_MAX_DIM(5,1),MIN_MAX_DIM(5,2)]*0.02)
    ylim([min(yrang), max(yrang)])
    set(gca, 'YDir', 'normal');  % Reverse the y-axis direction
    colorbar
    colormap(turbo)
    logticks = [1e-12, 1e-11, 1e-10,1e-9, 1e-8,1e-7, 1e-6,1e-5, 1e-4,1e-3, 1e-2, 1e-1];
    set(colorbar, 'Ticks', log10(logticks), 'Ticklabels', logticks)
    
    
    hold off
    title({"Z = 102, A = 256", "Energy = " + E(f) + " MeV", "Fusion chance: " + fusionChance(f) + "%"})
    xlabel("Mass assymetry α")
    ylabel("Quadropole deformation Q2")
end


function coord = coordFromLinearIdxShifted(idx, dimSize, MIN_MAX_DIM)
    tempidx = idx;
    multfactor = 1;
    
    for i = 1:length(dimSize) 
        multfactor = multfactor * dimSize(i);
    end 
    coord = zeros([size(dimSize),1]);
    for i = length(dimSize):-1:1
        if (i == 1)
            coord(i) = tempidx;
        else
            multfactor = floor(multfactor / dimSize(i));
            coord(i) = floor(round(tempidx - 1) / multfactor + 1);
            tempidx = floor(mod(tempidx - 1, multfactor) + 1);
        end
    end

    %for i = 1 : length(dimSize)
    %    coord(i) = coord(i) - 1 + MIN_MAX_DIM(i,1);
    %end
end 

function q2 = get_q2_val(AZ,AA,I_val)

    Aref = 240;
    Zref = 94;
    
    % Quadrupole deformation values for 240Pu - corresponding to index II.
    Q2vals = [4.81253,   6.82805,   8.90879,  11.06605,  13.31188,  15.65945,  18.12340,  20.72011,  23.46814,  26.38871,  29.50621,  32.84902,  36.45030,  38.39975,  40.34920,  42.47073,  44.59226,  46.91380,  49.23533,  51.79066,  54.34598,  57.17638,  60.00677,  63.12872,  66.25067,  69.48992,  72.72917,  76.04370,  79.35822,  86.12054,  92.99791,  99.97944, 107.05082, 114.19157, 121.42655, 128.69212, 136.03868, 143.41393, 150.81228, 158.20718, 165.78249, 173.09174, 180.81321, 188.18743, 195.65203, 203.06961, 210.66009, 218.12864, 225.62603, 233.36657, 240.74626, 247.79498, 254.57527];

    Q2factor = (AZ/Zref*AA/Aref).^(2/3);
    
    q2_scale = 0.01*3*AZ*1.2^2* AA^(2/3)/(4.0*pi) ;

    q2 = Q2factor*Q2vals(I_val)/q2_scale;
    
end