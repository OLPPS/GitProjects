function [GateDataFile_plusCS, CScase_plusCS]  = DetermineInducedPseudoHits(thresholdx, thresholdy, GateDataFile, xpixIDs, ypixIDs, xpos, ypos, xpixnum, ypixnum)
PixelPitch = 1; %This line allows the function to work in normalised mode, without having to specify that on the input, whilst still keeping the ocde clearly and easily generalisable for future developments.
%% Identify events to be assessed for charge sharing
%Set an array to track number of charge sharing events
CScase = uint16(zeros(size(GateDataFile,1),1));
%For events near left edge add 100
CScase(GateDataFile(:,xpos)<=thresholdx) = CScase(GateDataFile(:,xpos)<=thresholdx) + 100;
%For events near right edge add 1000
CScase((PixelPitch-GateDataFile(:,xpos))<=thresholdx) = CScase((PixelPitch-GateDataFile(:,xpos))<=thresholdx) + 1000;
%For events near bottom edge add 1
CScase(GateDataFile(:,ypos)<=thresholdy) = CScase(GateDataFile(:,ypos)<=thresholdy) + 1;
%For events near top edge add 10
CScase((PixelPitch-GateDataFile(:,ypos))<=thresholdy) = CScase((PixelPitch-GateDataFile(:,ypos))<=thresholdy) + 10;



%% Determine new array size and the coordinates for existing events to be mapped to
%Estimate how many charge sharing events will be needed
NumberOfSharedEvents = sum(ismember(CScase,[1000, 100, 10, 1] )) + 3*sum(ismember(CScase,[ 1010, 1001, 0110, 0101 ]));

%logical masks for edges and corners
EdgeCodes = [1000, 100, 10, 1];
CornerCodes = [1010, 1001, 0110, 0101];
IsAnEdge = ismember(CScase,EdgeCodes);
IsACorner = ismember(CScase,CornerCodes);

%Set up array for how many extra events will be needed
ExtraPixelsNeededForCS = ones(size(CScase));
%Add 1 event at edges
ExtraPixelsNeededForCS(IsAnEdge) = 2;
%Add 3 event at corners
ExtraPixelsNeededForCS(IsACorner) = 4;

%Calculate where events should go in new list
NewEventIndexes = cumsum(ExtraPixelsNeededForCS) - ExtraPixelsNeededForCS + 1;


%% Construct new array to hold events and populate with data from existing events
%Create new array to store the new CS events
GateDataFile_plusCS = zeros(NewEventIndexes(end),7);

%Transfer events over
GateDataFile_plusCS(NewEventIndexes,1:5) = GateDataFile;
GateDataFile_plusCS(NewEventIndexes,xpixnum) = xpixIDs;
GateDataFile_plusCS(NewEventIndexes,ypixnum) = ypixIDs;


%Copy over relevant CScase codes to identify which neighbouring pixels will
%experience CS
CScase_plusCS(NewEventIndexes) = CScase(:);

clear CScase IsAnEdge IsACorner GateDataFile xpixIDs ypixIDs


%% Add in the information for created pixels based on their previous events (DO NOT SHIFT IN SPACE OR TIME)
%For edges
[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),EdgeCodes);
GateDataFile_plusCS(targetIndices,1:7) = GateDataFile_plusCS(sourceIndices,1:7);

%For corners
[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),CornerCodes);
GateDataFile_plusCS(targetIndices,1:7) = GateDataFile_plusCS(sourceIndices,1:7);
targetIndices = targetIndices+1;
GateDataFile_plusCS(targetIndices,1:7) = GateDataFile_plusCS(sourceIndices,1:7);
targetIndices = targetIndices+1;
GateDataFile_plusCS(targetIndices,1:7) = GateDataFile_plusCS(sourceIndices,1:7);

%% Set CS map codes for all events, based on pixel in prototype CS map

%Set all CS map locations to 5 (central pixel) for original events
CScase_plusCS(NewEventIndexes(ExtraPixelsNeededForCS==1)) = 5;

%Set all CS map locations
[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),EdgeCodes(3));
CScase_plusCS(targetIndices) = 2;
CScase_plusCS(sourceIndices) = 5;
%mask10 pixels

[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),EdgeCodes(2));
CScase_plusCS(targetIndices) = 4;
CScase_plusCS(sourceIndices) = 5;
%mask100 pixels

[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),EdgeCodes(1));
CScase_plusCS(targetIndices) = 6;
CScase_plusCS(sourceIndices) = 5;
%mask1000 pixels

[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),EdgeCodes(4));
CScase_plusCS(targetIndices) = 8;
CScase_plusCS(sourceIndices) = 5;
%mask1 pixels


%Set all CS map locations for corners
[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),CornerCodes(3));
CScase_plusCS(targetIndices) = 1;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 2;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 4;
CScase_plusCS(sourceIndices) = 5;
%mask110 pixels

[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),CornerCodes(1));
CScase_plusCS(targetIndices) = 3;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 6;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 2;
CScase_plusCS(sourceIndices) = 5;
%mask1010 pixels

[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),CornerCodes(4));
CScase_plusCS(targetIndices) = 7;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 4;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 8;
CScase_plusCS(sourceIndices) = 5;
%mask101 pixels

[sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(CScase_plusCS(:),CornerCodes(2));
CScase_plusCS(targetIndices) = 9;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 6;
targetIndices = targetIndices + 1;
CScase_plusCS(targetIndices) = 8;
CScase_plusCS(sourceIndices) = 5;
%mask1001 pixels

%%Adjust PIXNUMS to account for new location of CS events
GateDataFile_plusCS(CScase_plusCS(:) == 1, xpixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 1, xpixnum) - 1;
GateDataFile_plusCS(CScase_plusCS(:) == 1, ypixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 1, ypixnum) + 1;

GateDataFile_plusCS(CScase_plusCS(:) == 2, ypixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 2, ypixnum) + 1;

GateDataFile_plusCS(CScase_plusCS(:) == 3, xpixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 3, xpixnum) + 1;
GateDataFile_plusCS(CScase_plusCS(:) == 3, ypixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 3, ypixnum) + 1;

GateDataFile_plusCS(CScase_plusCS(:) == 4, xpixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 4, xpixnum) - 1;

GateDataFile_plusCS(CScase_plusCS(:) == 6, xpixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 6, xpixnum) + 1;

GateDataFile_plusCS(CScase_plusCS(:) == 7, xpixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 7, xpixnum) - 1;
GateDataFile_plusCS(CScase_plusCS(:) == 7, ypixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 7, ypixnum) - 1;

GateDataFile_plusCS(CScase_plusCS(:) == 8, ypixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 8, ypixnum) - 1;

GateDataFile_plusCS(CScase_plusCS(:) == 9, xpixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 9, xpixnum) + 1;
GateDataFile_plusCS(CScase_plusCS(:) == 9, ypixnum) = GateDataFile_plusCS(CScase_plusCS(:) == 9, ypixnum) - 1;

%CScase_plusCS = CScase_plusCS'; %Convert from a row vector to a column vector
end


