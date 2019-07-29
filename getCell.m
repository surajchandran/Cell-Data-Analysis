function [pointsmax, optimalCellCount, optimalPmax, optimalMassDelta] = getCell(n,cells,fit)
    %Reading in cell data from cellData.xlsx and surface fits to Sam's
    %analysis from sweepfit.mat, if not loaded already. I fit points, brake time, endurance energy output, and
    %endurance time to the independent variables of accumulator mass delta and
    %max power usage.
    if nargin == 1
        celldata = xlsread('cellData.xlsx');
        sweepfit = load('sweepfit');
    else    
        celldata = cells;
        sweepfit = fit;
    end
    
    %Extracting the imported fit coefficients
    a = sweepfit.pointsfit(1);
    b = sweepfit.pointsfit(2);
    c = sweepfit.pointsfit(3);
    
    j = sweepfit.braketimefit(1);
    k = sweepfit.braketimefit(1);
    l = sweepfit.braketimefit(1);
    
    q = sweepfit.energyfit(1);
    r = sweepfit.energyfit(2);
    s = sweepfit.energyfit(3);
    
    x = sweepfit.timefit(1);
    y = sweepfit.timefit(2);
    z = sweepfit.timefit(3);
    
    %Reading the data of cell number n in the spreadsheet
    celln = celldata(n,:);
    
    %Extracting relevant cell data parameters
    cellV = (celln(2) + celln(3))/2;
    cellir = celln(4)/1000;
    cellE = celln(5);
    cellMass = celln(6)/1000;
    cellMaxDischarge = celln(11);
    cellMaxCharge = celln(13);
    
    %Setting independent variables of maxPower and massdelta
    maxPower = 10000:100:80000;
    massdelta = -25:0.5:25;
    
    %Setting up matrices for each parameter: endurance energy output,
    %endurance time, total dynamic points, RMS Power, brake time, and pack
    %mass. Matrices are m x n, where m is the number of MaxPower terms and
    %n is the number of massdelta terms to vary over.
    powersize = size(maxPower);
    masssize = size(massdelta);
    EndEOut = zeros(powersize(2),masssize(2));
    t = zeros(powersize(2),masssize(2));
    totalDynamicPts = zeros(powersize(2),masssize(2));
    Prms = zeros(powersize(2),masssize(2));
    braketime = zeros(powersize(2),masssize(2));
    packmass = zeros(powersize(2),masssize(2));
    
    %Using the surface fits obtained from Sam's lapsim data to calculate
    %each parameter value for each MaxPower/massdelta pair.
    for pindex = 1:powersize(2)
        for mindex = 1:masssize(2)
            EndEOut(pindex,mindex) = q*massdelta(mindex)+r*sqrt(maxPower(pindex))+s;
            t(pindex,mindex) = x*massdelta(mindex)+y*maxPower(pindex)^(-0.5)+z;
            totalDynamicPts(pindex,mindex) = a*massdelta(mindex)+b*sqrt(maxPower(pindex))+c;
            Prms(pindex,mindex) = -0.0000023*maxPower(pindex)^2 + 0.479*maxPower(pindex)+1949;
            braketime(pindex,mindex) = j*massdelta(mindex)+k*sqrt(maxPower(pindex))+l;
            packmass(pindex,mindex) = 26.8+massdelta(mindex); %26.8 kg is the current year's pack total cell mass
        end
    end    
    
    %Determining matrices for total number of cells and total energy
    %required. 
    cellCount = packmass/cellMass;
    %Determining total energy required. Formula used is endurance output
    %energy + internal resistive losses - regeneration.
    Etotal = EndEOut+(cellir/(cellV^2))*(cellCount.^(-1)).*(Prms.^2).*t - cellMaxCharge*cellV*cellCount.*braketime;
    
    %Setting up a matrix to hold points values for each MaxPower/massdelta
    %pair
    points = zeros(powersize(2),masssize(2));
    seriescells = 144;
    for pindex = 1:powersize(2)
        for mindex = 1:masssize(2) 
            if mod(cellCount(pindex,mindex),seriescells) < 1 || mod(cellCount(pindex,mindex),seriescells)>(seriescells-1)
                %Testing for a pack configuration at a particular voltage.
                %This can be adjusted by changing the value of the variable
                %"seriescells" in line 82.
                if cellCount(pindex,mindex)*3600*cellE > Etotal(pindex,mindex) && cellMaxDischarge*cellV*cellCount(pindex,mindex) > Prms(mindex)
                    %Testing to see if a particular MaxPower/massdelta pair has the
                    %capability to hit the Etotal calculated for that pair, and if
                    %it can continuously provide the Prms needed. If either test fails, the
                    %points for that run is set to 0. If both are true, then points
                    %is equal to the corresponding value in the previously
                    %calculated matrix.
                    points(pindex,mindex) = totalDynamicPts(pindex,mindex);
                else
                    points(pindex,mindex) = 0;
                end
            else
                points(pindex,mindex) = 0;
            end
        end
    end
    %Return the overall points maximum as a function output
    pointsmax = max(points(:));
    
    %Find the matrix indices of the max points value and extract the other
    %parameters with the same indices.
    if pointsmax > 0
        maxindex = find(points == pointsmax);
        maxindrow = mod(maxindex,powersize(2));
        if maxindrow == 0
            maxindrow = powersize(2);
        end
        maxindcolumn = ceil(maxindex/powersize(2));
        optimalCellCount = cellCount(maxindrow, maxindcolumn);
        optimalPmax = maxPower(maxindrow);
        optimalMassDelta = massdelta(maxindcolumn);
    else
        optimalCellCount = 0;
        optimalPmax = 0;
        optimalMassDelta = 0;
    end
end