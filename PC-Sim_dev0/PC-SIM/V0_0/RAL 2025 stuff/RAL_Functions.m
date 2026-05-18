close all
clc
clear


%%First load in the data

[DataFile, NumberOfEnergyBins, PixelsInX, PixelsInY, ScansPerFile, ScanFiles] = LoadData('Variables_HF15keV.mat');

for iijj = 1:1:ScanFiles
    FileFrame = squeeze(DataFile.data(:,:,:,:,iijj));
    for ij = 1:1:ScansPerFile
        FigNum = 0;

        SpectralFrame = squeeze(FileFrame(:,:,:,ij));
        EnergyFloor = 250
        EnergyCeiling = 600

        if EnergyCeiling <= EnergyFloor
            EnergyCeiling = NumberOfEnergyBins;
        end
        
%%Second, show all counts and all energy spectra
FigNum = AnalyseFrame(SpectralFrame, EnergyFloor, EnergyCeiling, FigNum, 'south');

%%Locating beam, isolate area around beam which it will sweep, with 2
%%pixels spare either side.
SweepInY = -5;
SweepInX = -5;
Border = 2;

if iijj == 1 && ij == 1
    [BoxedFrame, LowerY, UpperY, LowerX, UpperX] = BoxAroundMaxCounts(SpectralFrame,EnergyFloor, EnergyCeiling, SweepInX, SweepInY, Border);
else
    BoxedFrame = SpectralFrame(:,LowerY:UpperY,LowerX:UpperX);
end

FigNum = AnalyseFrame_Thresholded(BoxedFrame, EnergyFloor, 1024, FigNum, 'north');

BoxedFrame_stored(:,:,:,ij,iijj) = BoxedFrame;

    end
end

%%

for iijj = 1:10:ScanFiles
    for ij = 1:10:ScansPerFile
        FigNum = 0;
        dataFrame = squeeze(BoxedFrame_stored(:,:,:,ij,iijj));
FigNum = AnalyseFrame(dataFrame, EnergyFloor, EnergyCeiling, FigNum, 'north');
    end
end

%%
load("BoxedFrame_stored_HF15keV.mat")
%%
EnergyFloor = 250;
EnergyCeiling = 450;

sz = size(BoxedFrame_stored);  % [1024, 10, 10, 51, 51]
numProj = sz(4) * sz(5);  % Total number of projections (2601)

% Reshape A to combine last two dimensions
BoxedFrame_stored_timeOrdered = reshape(BoxedFrame_stored, [sz(1), sz(2), sz(3), numProj]);  % (1024,10,10,2601)

BoxedFrame_stored_timeOrdered_energySummed = squeeze(sum(BoxedFrame_stored_timeOrdered(EnergyFloor:EnergyCeiling,:,:,:),1));

%%
% Plot mesh figures for each projection
for p = 1:numProj
    figure(1)
    surf(BoxedFrame_stored_timeOrdered_energySummed(:,:,p));
    title(['Projection ' num2str(p)]);
    xlabel('X');
    ylabel('Y');
    zlabel('Summed Value');
    view(0,90)
    pause(2)
end





%%

maxCnts = max(max(max(BoxedFrame_stored_timeOrdered(EnergyFloor:EnergyCeiling,:,:))))

for n = 1:2601
n
    SpectraInFrame = BoxedFrame_stored_timeOrdered(EnergyFloor:EnergyCeiling,:,:,n);

% Now plot 10x10 subplots showing histograms of each (10,10) position
for y = 3:8
    for x = 3:8
figure(2);
        subplot(6,6,(y-3)*6 + x-2);
        plot(SpectraInFrame(:,y,x));
        ylim([0,maxCnts/4])
        title(['(' num2str(y) ',' num2str(x) ')'], 'FontSize', 6);
        set(gca, 'XTick', [], 'YTick', []);
    end
end
sgtitle(['Histograms of 1024 Channels Across 10x10 Grid ' num2str(n)]);
%pause(1)

end


%% STACK BASED ON LOCATION
initialX = 8;
initialY = 8;

xOffset = 7;
yOffset = 7;

xID = initialX;
yID = initialY;

IntraPixelLocations = uint32(zeros(1024,10,10));

IntraPixelLocations_3x3 = uint32(zeros(3,3,1024,10,10));

for ScanNumber = 1:50
ScanNumber
    for ScanPosition = 1:50
        % Work out which pixel value to look at, xID
xID = initialX - (floor((ScanPosition+xOffset)/10));
        % Work out x position within pixel, ix
ix = rem(ScanPosition,10)+1;

% Work out which pixel value to look at, yID
yID = initialY - (floor((ScanNumber+yOffset)/10));

% Work out y position within pixel, iy
iy = rem(ScanNumber,10)+1;

        IntraPixelLocations(:,iy,ix) = IntraPixelLocations(:,iy,ix) + BoxedFrame_stored(:,yID,xID,ScanPosition,ScanNumber);
     figure(4)
   
    for dy = -1:1
    for dx = -1:1

        relativeY = 2+dy;
        relativeX = 2+dx;

        IntraPixelLocations_3x3(relativeY,relativeX,:,iy,ix) = squeeze(IntraPixelLocations_3x3(relativeY,relativeX,:,iy,ix)) + squeeze(BoxedFrame_stored(:,yID+dy,xID+dx,ScanPosition,ScanNumber));

        subplot(3,3,(dy+2-1)*3 + dx+2);
        plot(BoxedFrame_stored(250:EnergyCeiling,yID+dy,xID+dx,ScanPosition,ScanNumber))
        ylim([0 10])
        title(['(' num2str(y) ',' num2str(x) ')'], 'FontSize', 6);
        set(gca, 'XTick', [], 'YTick', []);
    end
end
sgtitle(['Histograms of 1024 Channels Across 10x10 Grid ' num2str(ScanNumber) ' ' num2str(ScanPosition)]);

end
    
    end


    %%

    sumsIntraPixelLocations = squeeze(sum(IntraPixelLocations(EnergyFloor:EnergyCeiling,:,:),1));
    figure(7)
    surf(sumsIntraPixelLocations)
%%
InterPixSpec = uint32(zeros(1024,10,10));
pos = 1;
scan = 1;

                figure(5)

for scan = 1:10
    for pos = 1:10
        for dscan = 0:0
            for dpos = 0:4
                InterPixSpec(:,11-pos,11-scan) = InterPixSpec(:,11-pos,11-scan) + BoxedFrame_stored(:,initialX-dpos,initialY-dscan,pos+dpos*10,scan+dscan*10);
                sumsInterPixSpec = squeeze(sum(InterPixSpec(EnergyFloor:EnergyCeiling,:,:),1));
                surf(sumsInterPixSpec)
                view(0,90)
                %pause(0.5)
            end
        end
    end
end

sumsInterPixSpec = squeeze(sum(InterPixSpec,1));
lineSumsInterPixSpec = squeeze(sum(sumsInterPixSpec,2));

figure(6)
plot(lineSumsInterPixSpec)


%%

function [DataFile, numberOfEnergyBins, pixelsInX, pixelsInY, scansPerFile,scanFiles] = LoadData(VariablesFile)

load(VariablesFile);
DataFile = matfile(infile, 'Writable', false);
fprintf('Data contains:\n   %i energy bins\n   %i pixels in X\n   %i pixels in Y\nDerived from\n   %i scans from each of %i files\n\n', numberOfEnergyBins, pixelsInX, pixelsInY, scansPerFile,scanFiles)

end



function [figNum] = AnalyseFrame(data, energyFloor, energyCeiling, figNum,northSouth)
allCountsFrame = squeeze(sum(data(energyFloor:energyCeiling,:,:),1));
figNum = showHeatMap(allCountsFrame, energyFloor, energyCeiling, figNum,northSouth);

figNum = HistogramCounts(allCountsFrame, figNum,northSouth)

allPixelsSpectra = squeeze(sum(data(energyFloor:energyCeiling,:,:), [3, 2]));
figNum = showRawSpectrum(allPixelsSpectra,energyFloor,energyCeiling,figNum,northSouth);

end



function [figNum] = AnalyseFrame_Thresholded(data, energyFloor, energyCeiling, figNum,northSouth)
allCountsFrame = squeeze(sum(data(energyFloor:energyCeiling,:,:),1));
figNum = showHeatMap(allCountsFrame, energyFloor, energyCeiling, figNum,northSouth);

figNum = HistogramCounts(allCountsFrame, figNum,northSouth)

% 1. Extract the subset of the data within your energy range
dataSubset = data(energyFloor:energyCeiling, :, :);  % [energy, y, x]

% 2. Reshape to [energy, numPixels]
data2D = reshape(dataSubset, size(dataSubset,1), []);  % [energyBins, y*x]

% 3. Sum over energy for each pixel
pixelSums = sum(data2D, 1);  % [1, numPixels]

% 4. Identify pixels above threshold
validIdx = pixelSums > 10;  % [1, numPixels]

% 5. Keep only those pixels
validSpectra = data2D(:, validIdx);  % [energyBins, numValidPixels]

% 6. Sum over valid pixels
ThresholdedPixelsSpectra = sum(validSpectra, 2);  % [energyBins x 1]
figNum = showRawSpectrum(ThresholdedPixelsSpectra,energyFloor,energyCeiling,figNum,northSouth);

end






function [figNum] = showHeatMap(data, energyFloor, energyCeiling, figNum,northSouth)
figure(figNum+1)
figNum = figNum+1;
imagesc(data)
view(0,90)
guiPos = strcat(northSouth,'west');
movegui(guiPos)

colorbar;

% Label the axes
xlabel('X pixel');
ylabel('Y pixel');
title('Heatmap of pixel counts');

% Set axis properties
axis equal;
hold on;

% Set up grid
[numRows, numCols] = size(data);

% Add grid lines
for i = 1:numRows
    line([0.5, numCols+0.5], [i+0.5, i+0.5], 'Color', 'k'); % Horizontal lines
end
for j = 1:numCols
    line([j+0.5, j+0.5], [0.5, numRows+0.5], 'Color', 'k'); % Vertical lines
end

hold off;
end



function [figNum] = HistogramCounts(allCountsFrame, figNum,northSouth)
figure(figNum+1)
figNum = figNum+1;
histogram(allCountsFrame(:),'NumBins',100)
movegui(northSouth)
xlabel('Counts')
ylabel('Number of pixels')
title('Histogram of pixel counts (100 bins)')
end



function [figNum] = showRawSpectrum(data,energyFloor, energyCeiling, figNum,northSouth)
figure(figNum+1)
figNum = figNum+1;
plot([energyFloor:energyCeiling],data,'-r')
hold on
plot([energyFloor:energyCeiling],data,'ok')
hold off
xlim([energyFloor energyCeiling]);
guiPos = strcat(northSouth,'east');
movegui(guiPos)
xlabel('Channel number')
ylabel('Counts')
end


function [boxedData,lowerY,upperY,lowerX,upperX] = BoxAroundMaxCounts(data,energyFloor, energyCeiling, sweepInX, sweepInY, border)


ThresholdedCountsFrame = squeeze(sum(data(energyFloor:energyCeiling,:,:),1));

% Find maximum value and its linear index
[maxValue, linearIndex] = max(ThresholdedCountsFrame(:));

% Convert linear index to row and column indices
[row, col] = ind2sub(size(ThresholdedCountsFrame), linearIndex);

if sweepInY < 1
    row = row + sweepInY;
    sweepInY = - sweepInY;
end

if sweepInX < 1
    col = col + sweepInX;
    sweepInX = - sweepInX;
end

lowerY = row-border;
upperY = row+sweepInY+border;
lowerX = col-border;
upperX = col+sweepInX+border;

while lowerY < 1
    lowerY = lowerY+1;
    upperY = upperY+1;
end
while upperY > size(ThresholdedCountsFrame,1)
    upperY = upperY-1;
    lowerY = lowerY-1;
end

while lowerX < 1
    lowerX = lowerX+1;
    upperX = upperX+1;
end
while upperX > size(ThresholdedCountsFrame,2)
    upperX = upperX-1;
    lowerX = lowerX-1;
end

boxedData = data(:,lowerY:upperY,lowerX:upperX);

end