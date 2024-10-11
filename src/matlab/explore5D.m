clear
close
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

MIN_MAX_DIM = [MIN_II, MAX_II; MIN_JJ, MAX_JJ; MIN_KK, MAX_KK; MIN_LL, MAX_LL; 0, MAX_MM];
dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ, 1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM];
%MIN_MAX_DIM = [MIN_II, MAX_II; MIN_JJ, MAX_JJ; MIN_KK, MAX_KK; MIN_LL, MAX_LL; MIN_MM, MAX_MM];
%dimSize = [1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ, 1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM - MIN_MM];
pd = readmatrix("../../data/PD", "Delimiter", ";"); %raw prob distribution
startCoord = [17,7,3,13,30];
%pd = str2double(pdstr)
%This is very wrong.

A = zeros(dimSize);

for i = 1:length(pd)
    if(mod(i,100000) == 0)
        disp(100*i/length(pd) + " percent complete")
    end
    c = coordFromLinearIdxShifted(i, dimSize, MIN_MAX_DIM);
    A(c(1),c(2),c(3),c(4),c(5)) = pd(i);
end
%%
close all
figure
TEMP = (A);
TEMP(find(TEMP<1e-20)) = 1e-20;
TEMP = log10(TEMP);
TEMP = permute(TEMP,[1,5,2,3,4]);

for i = 1:dimSize(2)
    imagesc(TEMP(:,:,i,3,13));
    caxis([min(TEMP,[],'all'), max(TEMP,[],'all')]);
    colorbar
    
    title("Neck: coordinate: " + num2str(i))
    xlabel("Mass assymetry")
    ylabel("Elongation Q")
    pause(1)
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