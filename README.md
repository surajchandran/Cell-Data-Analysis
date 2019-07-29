# Cell-Data-Analysis

The Matlab data file "sweepfit" has interpolation information for the massdelta/endurancepower sweep run by Sam.
getCell optimizes points for a single cell. analyzeCellData runs getCell across all cells in the Cell Data excel file.
The output of analyzeCellData is a matrix array.
Column 1: cell index in the excel file.
Column 2: Optimized points.
Column 3: optimal cell count.
Column 4: optimal max endurance power.
Column 5: optimal mass delta.
After pulling the files in the repository, simply run analyzeCellData in the Matlab console.
