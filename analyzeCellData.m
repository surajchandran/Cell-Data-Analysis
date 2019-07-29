function pointstable = analyzeCellData()
    %Reading cell data and loading surface fits to Sam's lapsim analysis
    %with independent variables of MaxPower and massdelta.
    cells = xlsread('cellData.xlsx');
    fit = load('sweepfit');
    
    %Setting up a points table. Tabulates points for each cell in the
    %format [vaultindex score].
    pointstable = zeros(81,5);
    for i=1:81
       [points, optCellCount, optPower, optimalMass] = getCell(i,cells,fit);
       pointstable(i,:) = [i points optCellCount, optPower, optimalMass];
    end
end