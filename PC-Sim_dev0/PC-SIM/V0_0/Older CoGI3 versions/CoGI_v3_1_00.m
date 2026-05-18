%% Function which performs the main bulk of CoGI
% ParamsFile holds various parameters needed at this stage, e.g. which CSCAs to use.
% GeometryFile defines the geometry of the detector with respect to the x-ray source in STANDARDISED format.
% NoiseParametersFile imports information regarding the noise profile of the detector.
% MCDataMatFileName contains the MC data in STANDARDISED format.
% t1, and t2 are the start and end times of the slice being processed.
% RANDOMSEED is set to allow reproducibillity of the simulation and fair comparisson between various CSCAs or ACTS etc.
% CIEtolerance is a threshold which can be set such that events with a CIE value below this point are not processed further.
% SectionsToRun this is a vector which stores which sections of CoGI to run. See documentation for more details.
% MapAddresses is an array which stores the addresses of the CIE maps for centre, edge, corner, CSEdge and CSCorner in STANDARISED format, in that order.
% VERBOSE sets text output (1 = a lot, default = basic)
% DEBUGGING sets whether debugging mode is on (see documentation for more details.)
% Label is the base file name for outputs.
% THRESHOLDS is a vector containing the thresholds to be applied in keV.

function CoGI_v3_1_0(ParamsFile,GeometryFile,NoiseParametersFile,MCDataMatFileName,t1,t2, RANDOMSEED,CIEtolerance, SectionsToRun, MapAddresses, SigmaTolerance, VERBOSE,DEBUGGING,Label,THRESHOLDS)

NumOfThresh = length(THRESHOLDS);

if VERBOSE == 1
startTime = datetime
end
if SectionsToRun(1) == 1
CompressArray = 1;
adjustForDriftTime = 0;

%Initialise debugging/error flags
%[Flags, IndicesChecker] = InitialiseDebuggingFlags(); %This line could be removed?

%% Precompute all CIE maps needed and save to files
%{
    MapLabels = {'Central','Edge','Corner','CSLeftPixel','CSTopLeftPixel'};

    for MapPos = 1:5
        ComsolFileAddress = strcat(COMSOLBaseFileName,MapLabels{MapPos});
        [CIEMap, CIETimeChosenIndices,xVals, yVals, zVals] = GenerateCIEMap(ComsolFileAddress,'Max');
        fullFilePath = fullfile('CIEMaps_', ComsolFileAddress);
        save(fullFilePath, 'CIEMap', 'CIETimeChosenIndices','xVals', 'yVals', 'zVals', 'ComsolFileAddress');
        MapAddresses{MapPos} = fullFilePath;
    end
    SelectedSignalStructured3DMap
end
if VERBOSE == 1
CIEMaps_Computed = datetime
%}
%% Load Parameters for detector and simulation
%%Load processing parameters
load(ParamsFile)

%%Load simulation geometry
load(GeometryFile)

%%Set geometry parameters with easier to read variable names
            PixelWidth = GEOMETRIES(1); %mm
            PixelPitch = GEOMETRIES(2); %mm
            PixelDepth = GEOMETRIES(3); %mm
            xOffset = GEOMETRIES(4); %mm x translation of detector centre from origin in x
            yOffset = GEOMETRIES(5); %mm y translation of detector from origin in y
            zOffset = GEOMETRIES(6); %mm z translation of detector from origin in z
            MagXEdgeSpace = GEOMETRIES(7); %mm Shifts to correct for symmetric expansion of detector in MC, ignoring translations. Array width in x /2 if detector centred on origin
            MagYEdgeSpace = GEOMETRIES(8); %mm Shifts to correct for symmetric expansion of detector in MC, ignoring translations. Array length in y /2 if detector centred on origin
            PixelsInXDirection = GEOMETRIES(9); %pixels - 1
            PixelsInYDirection = GEOMETRIES(10); %pixels - 1

%%Calculate Noise maps
NoisinessMaps_FileName = BuildDetector(NoiseParametersFile,GEOMETRIES(9),GEOMETRIES(10),RANDOMSEED);
load(NoisinessMaps_FileName)


if VERBOSE == 1
NoiseMaps_Computed = datetime
end


%% Process MC data based on Geometry

GateDataFile = ExtractMCDataTimeSlice(MCDataMatFileName, t1, t2);

%%Set out column labels
            Columns = 11;
            xpos = 1;
            ypos = 2;
            zpos = 3;
            time = 4;
            energy = 5;
            xpixnum = 6;
            ypixnum = 7;
            pixtype = 8;
            CIE = 9;
            Signal = 10;
            DriftState = 11;

%%Identify zero energy or out of range events
ZeroEnergyMask = GateDataFile(:,energy)==0;
energyzeroes = sum(ZeroEnergyMask);

MaxUncorrectedX = (PixelPitch * (PixelsInXDirection))+ xOffset - MagXEdgeSpace;
MaxUncorrectedY = (PixelPitch * (PixelsInYDirection))+ yOffset - MagYEdgeSpace;
MinUncorredtedX = 0 + xOffset - MagXEdgeSpace;
MinUncorredtedY = 0 + yOffset - MagYEdgeSpace;

ValidPixels = ( (GateDataFile(:,xpos) >= MinUncorredtedX) & (GateDataFile(:,ypos) >= MinUncorredtedY) & (GateDataFile(:,xpos) <= MaxUncorrectedX) & (GateDataFile(:,ypos) <= MaxUncorrectedY) );
ValidPixels = ValidPixels & ~ZeroEnergyMask;
OutsideStrikes = sum(~ValidPixels);
clear MaxUncorrectedX MaxUncorrectedY MinUncorrectedX MinUncorrectedY ZeroEnergyMask


%%Adjust Global Coordinates for offsets and origin centred nature of MC coordinates...
deltaX = MagXEdgeSpace- xOffset; %MagDEdgeSpace accounts for multiple pixels in this direction. For the 1 pixel case it reduces to pixel pitch/2
deltaY = MagYEdgeSpace - yOffset;
deltaZ = PixelDepth/2.0 - zOffset;
GateDataFile(ValidPixels,xpos) = GateDataFile(ValidPixels,xpos) + deltaX;
GateDataFile(ValidPixels,ypos) = GateDataFile(ValidPixels,ypos) + deltaY;
            if SourceDetectorOrientation == 'rev'
                GateDataFile(ValidPixels,zpos) = PixelDepth - (GateDataFile(ValidPixels,zpos) + deltaZ); %This is to account for z being reversed between Gate and Comsol in 'rev' orientation
            elseif SourceDetectorOrientation == 'std'
                GateDataFile(ValidPixels,zpos) = GateDataFile(ValidPixels,zpos) + deltaZ;
            else 
                fprintf('Unrecognised Source Detector Orientation Specified.')
                quit;
            end

%%Determine pixel numbers etc.
xpixIDs = zeros(size(GateDataFile,1),1);
ypixIDs = zeros(size(GateDataFile,1),1);

xpixIDs(ValidPixels) = fix(GateDataFile(ValidPixels,xpos) / PixelPitch)+1;
ypixIDs(ValidPixels) = fix(GateDataFile(ValidPixels,ypos) / PixelPitch)+1;

%%Calculate intrapixel coordinates
GateDataFile(ValidPixels,xpos) = GateDataFile(ValidPixels,xpos) - (xpixIDs(ValidPixels)-1)*PixelPitch; %-1 to correct for the nth pixel only needing to cut (n-1) previous pixels
GateDataFile(ValidPixels,ypos) = GateDataFile(ValidPixels,ypos) - (ypixIDs(ValidPixels)-1)*PixelPitch; %-1 to correct for the nth pixel only needing to cut (n-1) previous pixels

%%Normalise intrapixel coordinates (needed to avoid CS map induced location mismatches. See Note 1).
GateDataFile(ValidPixels,xpos) = GateDataFile(ValidPixels,xpos) / PixelPitch;
GateDataFile(ValidPixels,ypos) = GateDataFile(ValidPixels,ypos) / PixelPitch;
GateDataFile(ValidPixels,zpos) = GateDataFile(ValidPixels,zpos) / PixelDepth;


if VERBOSE == 1
CoreMCDataProcessing_Computed = datetime - NoiseMaps_Computed
end


%% Implement CS module

%Set which parts of the map should be considered for charge sharing. The
%default for point-charge clouds is the pixel street only. If finite Gaussian
%charge clouds are used, this number should be changed based on the drop in
%CIE in the CIE maps post convolution with the largest charge cloud size
%considered.
thresholdx = (PixelPitch - PixelWidth)/2.0;
thresholdy = (PixelPitch - PixelWidth)/2.0;

%%Normalise spatial threshold coordinates (needed to avoid CS map induced location mismatches. See Note 1).
thresholdx = thresholdx/PixelPitch;
thresholdy = thresholdy/PixelPitch;


[GateDataFile_plusCS, CScase_plusCS] = DetermineInducedPseudoHits(thresholdx,thresholdy,GateDataFile(ValidPixels,:),xpixIDs(ValidPixels),ypixIDs(ValidPixels),xpos,ypos,xpixnum,ypixnum);
clear GateDataFile xpixIDs ypixIDs ValidPixels


%Only consider CS induced events within the physical pixel array.
ValidEvents = (GateDataFile_plusCS(:,xpixnum) >= 1 & GateDataFile_plusCS(:,xpixnum) <= PixelsInXDirection & GateDataFile_plusCS(:,ypixnum) >= 1 & GateDataFile_plusCS(:,ypixnum) <= PixelsInYDirection);

if VERBOSE == 1
AddtionOfCSEvents_Computed = datetime
end


%% Caclulate CIE values
%Initialise CIE values
CIE_Values = zeros(size(GateDataFile_plusCS,1),1);

%%Classify pixels as central, edge or corner type: 0 = central, 1 = xedge, 2 = yedge, 3 = corner
PixTypeIDs = zeros(size(GateDataFile_plusCS,1),1);
xEdges = GateDataFile_plusCS(:,xpixnum) == 1 | GateDataFile_plusCS(:,xpixnum) == PixelsInXDirection;
PixTypeIDs(xEdges) = 1; %Can do this as only xEdges checked so far
yEdges = GateDataFile_plusCS(:,ypixnum) == 1 | GateDataFile_plusCS(:,ypixnum) == PixelsInYDirection;
PixTypeIDs(yEdges) = PixTypeIDs(yEdges) + 2; %So pixels in yEdges and xedges get total score of 2.

%Step through for each map type: Central, XEdges, YEdges, Corner
for PixelLocationType = 0:3
    %Load in map data
    switch PixelLocationType
        case 0
            CIEMapAddress = MapAddresses{1};
        case 1
            CIEMapAddress = MapAddresses{2};
        case 2
            CIEMapAddress = MapAddresses{2};
        case 3
            CIEMapAddress = MapAddresses{3};
    end
    load(CIEMapAddress);
    %%Normalise CIE map coordinates (needed to avoid CS map induced location mismatches. See Note 1).
    xVals = (xVals - min(xVals)) / (max(xVals) - min(xVals)); 
    yVals = (yVals - min(yVals)) / (max(yVals) - min(yVals)); 
    zVals = (zVals - min(zVals)) / (max(zVals) - min(zVals)); 

    %Identify core pixels using this map
    PixelsBeingCalculated = ( (PixTypeIDs(:) == PixelLocationType) & ValidEvents(:) & (CScase_plusCS(:) == 5) );

    %Build interpolation grids
    [X, Y, Z] = ndgrid(xVals, yVals, zVals);

    %if size(GateDataFile_plusCS(PixelsBeingCalculated,xpos),1) ~= 0 && size(GateDataFile_plusCS(PixelsBeingCalculated,ypos),1) ~= 0 && size(GateDataFile_plusCS(PixelsBeingCalculated,zpos),1) ~= 0
    %Trilinearly interpolate CIE values, allowing extrapolation for
    %rounding errors near edges
    F = griddedInterpolant({xVals, yVals, zVals}, CIEMap, 'linear', 'linear');
    CIE_Values(PixelsBeingCalculated) = F(GateDataFile_plusCS(PixelsBeingCalculated, zpos), GateDataFile_plusCS(PixelsBeingCalculated, ypos), GateDataFile_plusCS(PixelsBeingCalculated, xpos)); %HardCorrection for COMSOL to GATE mismatch directions. Add to documentation and fix after japan.
    %end

    clear xqs yqs zqs PixelsBeingCalculated X Y Z

end


%Identify CS pixels and the appropriate map to use based on their CScase values
for CSMap = 1:9
    %Load in map data and rotate as needed
    switch CSMap
        case 1
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 0), [3 1 2]);
        case 2
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 0), [3 1 2]);
        case 3
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 1), [3 1 2]);
        case 4
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 3), [3 1 2]);
        case 5
            CIEMapAddress = MapAddresses{1};
            load(CIEMapAddress);
        case 6
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 1), [3 1 2]);
        case 7
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 3), [3 1 2]);
        case 8
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 2), [3 1 2]);
        case 9
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 2), [3 1 2]);
    end

    %Identify core pixels using this map
    PixelsBeingCalculated = ( ValidEvents(:) & (CScase_plusCS(:) == CSMap) );

    %%Normalise CIE map coordinates (needed to avoid CS map induced location mismatches. See Note 1).
    xVals = (xVals - min(xVals)) / (max(xVals) - min(xVals)); 
    yVals = (yVals - min(yVals)) / (max(yVals) - min(yVals)); 
    zVals = (zVals - min(zVals)) / (max(zVals) - min(zVals)); 

    
%%FOR FIRST TEST RUN
% Check size alignment
assert(size(CIEMap,1) == length(xVals), 'Mismatch: xVals must match size(CIEMap,1)');
assert(size(CIEMap,2) == length(yVals), 'Mismatch: yVals must match size(CIEMap,2)');
assert(size(CIEMap,3) == length(zVals), 'Mismatch: zVals must match size(CIEMap,3)');

disp('CIEMap axis alignment with xVals, yVals, zVals is correct.');
%%FOR FIRST TEST RUN

  %Build interpolation grids
    [X, Y, Z] = ndgrid(xVals, yVals, zVals);

    %if size(GateDataFile_plusCS(PixelsBeingCalculated,xpos),1) ~= 0 && size(GateDataFile_plusCS(PixelsBeingCalculated,ypos),1) ~= 0 && size(GateDataFile_plusCS(PixelsBeingCalculated,zpos),1) ~= 0
    %Trilinearly interpolate CIE values, allowing extrapolation for
    %rounding errors near edges
    F = griddedInterpolant({xVals, yVals, zVals}, CIEMap, 'linear', 'linear');
    CIE_Values(PixelsBeingCalculated) = F(GateDataFile_plusCS(PixelsBeingCalculated, zpos), GateDataFile_plusCS(PixelsBeingCalculated, ypos), GateDataFile_plusCS(PixelsBeingCalculated, xpos));%HardCorrection for COMSOL to GATE mismatch directions. Add to documentation and fix after japan.
    %end
clear xqs yqs zqs PixelsBeingCalculated X Y Z

end

EventsToTransfer = ValidEvents & (abs(CIE_Values) >= CIEtolerance);

if VERBOSE == 1
InterpolatedCIEValues_Computed = datetime
end


%% Adjust for drift time if chosen
        if adjustForDriftTime == 1
        GateDataFile_plusCS(:,time) = GateDataFile_plusCS(:,time) + GateDataFile_plusCS(:,PixelDepth*(zpos))/(mu*Vb); %Time + DriftTime %CHECK if it's zpos or 1-zpos. zpos used atm, assuming z increases towards cathode (std orientation?)
        %Check the 1-zpos. Should it be xpos? or 1-xpos? or zpos etc.? Make a document to help show orientation%(d*l)/(mu*Vb) but l = 1 as normalised spatial coordinates are used.
        end

        if VERBOSE == 1
SignalDriftTimeAdjustments_Computed = datetime
end

%% First pass adjust to group results in same pixel and ns to reduce data input to and size of pixelated array
if CompressArray == 1
    xmax = PixelsInXDirection;
ymax = PixelsInYDirection;

Signals = CIE_Values(:).*GateDataFile_plusCS(:,energy)*1000;

for X = 1:xmax
    for Y = 1:ymax
        % Find all rows in A where xpix == X and ypix == Y
        pixelMask = (GateDataFile_plusCS(:,xpixnum) == X) & (GateDataFile_plusCS(:,ypixnum) == Y) & EventsToTransfer;
        pixelEvents = find(pixelMask);

        % Loop through sorted events and merge duplicates
        n = length(pixelEvents);
        while n > 1
            currIdx = pixelEvents(n);
            prevIdx = pixelEvents(n-1);

            if GateDataFile_plusCS(currIdx, time) == GateDataFile_plusCS(prevIdx, time)
                Signals(prevIdx) = Signals(prevIdx) + Signals(currIdx);
                CIE_Values(currIdx) = -89865;
            end
            n = n - 1;
        end
    end
end

EventsToTransfer = EventsToTransfer & CIE_Values(:)~=-89865;

end



%%ADD NOISE HERE FOR HOT PIXELS HERE IF USING HOT PIXEL GENERATED NOISE
%%INSTANCES HERE. NOT CLEAR THAT THIS IS NEEDED AT ALL THOUGH.

%% Build PixelatedBoard and other tools - during later OPTIMISATION, this could be done with a cell array to prevent excessive RAM usage in cases where some pixels are significantly hotter than others. This would also allow hot pixel nosie to be added.

Sz = sum(EventsToTransfer,1);

%Change in symbol to be compatible with older versions of ACTS and CSCAs.
XMax = PixelsInXDirection;
YMax = PixelsInYDirection;

%Calculate Bin sizes needed
Counter = zeros(XMax,YMax);

Counter(:,:) = accumarray([GateDataFile_plusCS(EventsToTransfer,xpixnum), GateDataFile_plusCS(EventsToTransfer,ypixnum)], 1, [XMax, YMax]);

MaxBinDepth = max(Counter(:))+2;

    if NOISEON == 0
        PixelatedBoard = zeros(XMax,YMax,4,MaxBinDepth); 
    elseif NOISEON == 1
        PixelatedBoard = zeros(XMax,YMax,5,MaxBinDepth); %%NOISE changed to 5
    end



excessiveTime = 2.0*(max(GateDataFile_plusCS(EventsToTransfer,time))+CSCA_SearchingTime);
PixelatedBoard(1:XMax,1:YMax,1,2:MaxBinDepth) = excessiveTime;
Headers = 2 * ones(XMax,YMax,2); %This ensures first value has time = 0; May need to account for this if doing a pixelated N ~= MasterList N test later? Equally means there will be NumOfPixels more events than expected if using headers(2) to calcualte events without subtracting 1 from headers(2) first.
MasterList(1:Sz,1:4) = -1000;

%PixelatedBoard labels
Tyme = 1; ProtoSig = 2; HashNum = 3; PrevSignal = 4;

if NOISEON == 1
NoiseAtTyme = 5; %NOISE This could be moved to replace HashNum if not in debugging mode during later OPTIMISATION
%%INITIAL NOISE ADDED HERE. N.B THAT NOISE IS ALREADY IN keV, SO DOES NOT NEED TO BE MULTIPLIED BY 1000
rng(RANDOMSEED+18+uint32(100000*excessiveTime), 'twister');
PixelatedBoard(:,:,5,:) = NoiseLevelMap_Adjusted .* randn(PixelsInXDirection, PixelsInYDirection, MaxBinDepth);
PixelatedBoard(:,:,2,:) = -1000000; %This adds a very negative signal to empty events, to prevent noise creating counts where events do not actually exist, as this is accounted for in a different module later on. See Note 2.
end

%Header labels - For reference in notes, not used.
Pointer = 1; ListSize = 2;

%MasterList
NNum = 1; %This can be removed during later OPTIMISATION: not used.
Xval = 2; Yval = 3; Tag = 4;

%% Build PixelatedBoard and other tools

IndicesOfEventsToTransfer = find(EventsToTransfer);

    for evnt = 1:length(IndicesOfEventsToTransfer)
        N = IndicesOfEventsToTransfer(evnt);
        X = GateDataFile_plusCS(N,xpixnum);
        Y = GateDataFile_plusCS(N,ypixnum);
        PixelatedBoard(X,Y,1,Headers(X,Y,2)) = GateDataFile_plusCS(N,time);
        PixelatedBoard(X,Y,2,Headers(X,Y,2)) = Signals(N);
        PixelatedBoard(X,Y,3,Headers(X,Y,2)) = N; 
        %   PixelatedBoard(X,Y,4,Headers(X,Y,2)) = 1; %CAN BE REMOVED OPTIMISATION
        Headers(X,Y,2) = Headers(X,Y,2)+1;
    end


        %MasterList(:,1) = N; % Could be removed during OPTIMISATION? This would require a rewrite as shifts the index for MasterList.
        MasterList(:,2) = GateDataFile_plusCS(EventsToTransfer,xpixnum);
        MasterList(:,3) = GateDataFile_plusCS(EventsToTransfer,ypixnum);
        MasterList(:,4) = ones(Sz,1);

        Headers(:,:,2) = Headers(:,:,2) - 1;
        %MasterList_Backup = MasterList; %THIS CAN BE REMOVED DURING OPTIMISATION

%clear GateDataFile_plusCS

%%

if DEBUGGING == 1
    outputFileName = strcat(Label,'_DDD_PixBrdPopulated');
    save(outputFileName);
    varsBefore = who;
elseif DEBUGGING == 2
    outputFileName = strcat(Label,'_DDD_PixBrdPopulated_v7pt3');
    save(outputFileName, '-v7.3');
    varsBefore = who;
end

end

%Load data in if needed/indicated
if SectionsToRun(2) == 1
    if SectionsToRun(1) == 0
        if DEBUGGING == 1
            inputFileName = strcat(Label,'_DDD_PixBrdPopulated');
            load(Label_DDD_PixBrdPopulated)
        elseif DEBUGGING == 2
            inputFileName = strcat(Label,'_DDD_PixBrdPopulated_v7pt3');
            load(Label_DDD_PixBrdPopulated)
        end
    end


    %% Apply CSCAs and or ACTS as indicated, with or without noise.
    Thresholds = THRESHOLDS;

    for NoiseOff = 1:1
if NOISEON == 0

for NoiseOff_ApplyFDR = 1:1        
%% Apply FDR counting scheme if requested - In current form, this needs to be implemented before decay correction BUT ALSO IS NOT AFFECTED BY PULSE SHAPE IN THAT CASE
fdrstart = datetime;
for Applying_FDR_ACTs = 1
    if FDR == 1
        NoCSCA_FDR_Counters = zeros(XMax, YMax,NumOfThresh);
        NoCSCA_FDR_EventsProcessed = 0;
        NoCSCA_PercentageCompleted = 0;
        NoCSCA_FDR_updateInterval = Sz*0.05;
        NoCSCA_FDR_nextUpdate = NoCSCA_FDR_updateInterval;
        for X = 1: XMax
            for Y = 1:YMax

%%New FDR version starts here                
PrevSignal = 0;
pointer = 1;
FrameSituation = 0;
OutsideOfFrame = 0;
while pointer <= Headers(X,Y,2) %Known events for this pixel. While needed to allow continue to loop over events in the out of frame section
  
    %%Percentage complete tracker
    NoCSCA_FDR_EventsProcessed = NoCSCA_FDR_EventsProcessed + 1;
    if NoCSCA_FDR_EventsProcessed >= NoCSCA_FDR_nextUpdate
fprintf('FDR: %.2f%% completed\n', (NoCSCA_FDR_EventsProcessed/Sz)*100);
NoCSCA_FDR_nextUpdate = NoCSCA_FDR_nextUpdate + NoCSCA_FDR_updateInterval;
    end
    %%Percentage complete tracker

    if FrameSituation  == 0
        %DO REGULAR SIG CALC
        if OutsideOfFrame == 1
            deltaT = (PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1));
            PrevSignal = CurrentSignal*exp(-DecayConstant*deltaT); %calculate and record residual signal in pixel before current event
        end
        CurrentSignal = PixelatedBoard(X,Y,2,pointer) + PrevSignal; %add current event onto residual signal
        %COMP TO THRESHOLDS
        BelowThresh = 0;
        for Thresh = 1:NumOfThresh
            if PrevSignal <= Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                BelowThresh = Thresh;
                break
            end
        end
        if BelowThresh ~= 0
            AboveThresh = 0;
            for Thresh = BelowThresh:NumOfThresh
                if CurrentSignal > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                    AboveThresh = Thresh;
                else
                    break
                end
            end
            %IF THRESHOLD TRIGGERED
            if AboveThresh ~= 0
                %ADD TO COUNTERS
                NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                %CALC NEW TIMES
                ResetStartTime = PixelatedBoard(X,Y,1,pointer) + FDR_RESETTIME;
                ResetEndTime = ResetStartTime + FDR_RESETDURATION;
                %CHANGE FRAMESITUATION TO 1
                FrameSituation = 1;
            end
        end
        OutsideOfFrame = 1;
    else %if FrameSituation ~= 0, so In Frame
        %CALCULATE DELTAT
        deltaT = (PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1));
        %COMPARE DELTAT TO TIMES TO DETERMINE SITUATIONS
        %IF IN FRAME
        if PixelatedBoard(X,Y,1,pointer) < ResetStartTime
            PrevSignal = CurrentSignal*exp(-DecayConstant*deltaT); %calculate and record residual signal in pixel before current event
            %KEEP COUNTING AS USUAL
            %DO REGULAR SIG CALC
            CurrentSignal = PixelatedBoard(X,Y,2,pointer) + PrevSignal; %add current event onto residual signal

            %COMP TO THRESHOLDS
            BelowThresh = 0;
            for Thresh = 1:NumOfThresh
                if PrevSignal <= Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                    BelowThresh = Thresh;
                    break
                end
            end
            if BelowThresh ~= 0
                AboveThresh = 0;
                for Thresh = BelowThresh:NumOfThresh
                    if CurrentSignal > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                        AboveThresh = Thresh;
                    else
                        break
                    end
                end
                %IF THRESHOLD TRIGGERED
                if AboveThresh ~= 0
                    %ADD TO COUNTERS
                    NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                end
            end

            %ELESIF IN RESET
        elseif PixelatedBoard(X,Y,1,pointer) > ResetEndTime %(IMPLICIT BEYOND FRAME, also implicitly ignores events in reset period)
            %PREVSIGNAL = 0
            PrevSignal = 0;
            %FRAMESITUATION == 0
            FrameSituation = 0;
            OutsideOfFrame = 0;
            %BREAK
            continue
            %END
        end
    end
    pointer = pointer +1;
end
%%New FDR version ends
            end
        end
    end
end
fprintf('FDR: 100%% completed\n');

%THIS MODULE SHOULD BE TESTED IN A SANDBOX
  fdrstop = datetime;
  fdrduration = fdrstop - fdrstart

end

for NoiseOff_DecayCorrectSignals = 1:1
%% Decay correct signals to be a running total
switch non_impulsePulseShape
case 1
    
case 2
for X = 1: XMax
    for Y = 1:YMax
        for pointer = 2:Headers(X,Y,2)
            deltaT = PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1); %calculate time since last event
            PixelatedBoard(X,Y,4,pointer) = PixelatedBoard(X,Y,2,pointer-1)*exp(-DecayConstant*deltaT); %calculate and record residual signal in pixel before current event
            PixelatedBoard(X,Y,2,pointer) = PixelatedBoard(X,Y,2,pointer) + PixelatedBoard(X,Y,4,pointer); %add current event onto residual signal
        end
    end
end

if DEBUGGING == 1
fprintf('Saving Output Data For Debug Logs...');
startA = datetime;
save(strcat(LABEL,'_DBG4'), 'PixelatedBoard', '-v7.3')
stopB = datetime;
fprintf('Completed in (hr:min:sec): %s\n', stopB-startA);
end

case 3 
%%%%%%%%%%%%%%%%%%%%%%%%%
%3). From empirical data (slowest, empirical model)\n',[1,2,3],REPEATEDINCORRECTCHOICEERRORMESSAGE);
%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

A = datetime;
for NoiseOff_ApplyNonFDR_ACTS = 1:1
%% Apply ACTS
%Define Thresholds
algstart = datetime;
NumOfThresh = max(size(Thresholds));
% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
NoCSCA_STD_Counters = zeros(XMax, YMax,NumOfThresh);
if PCS
    NoCSCA_PCS_Counters = zeros(XMax, YMax,NumOfThresh);
end
if DLR
    NoCSCA_DLR_Counters = -1 * ones(XMax, YMax,NumOfThresh); %-1 as initialising the buffer causes a single count to go to each threshold, so this means the system activation only will return 0 counts
end
if SR
    NoCSCA_SR_Counters = -1 * ones(XMax, YMax,NumOfThresh); %-1 as initialising the buffer causes a single count to go to each threshold, so this means the system activation only will return 0 counts
end
if IDEAL
NoCSCA_IDEAL_Counters = zeros(XMax, YMax,NumOfThresh);
end

NoCSCA_EventsProcessed = 0;
for X = 1:XMax
    for Y = 1:YMax
        DLR_HighestThreshReached = NumOfThresh;
        DLR_closetime = - abs(2.0*PixelatedBoard(X,Y,1,2) + DLR_INTEGRATIONTIME);
        SR_closetime = - abs(2.0*PixelatedBoard(X,Y,1,2) + SR_INTEGRATIONTIME);
        SR_Buffer(1:NumOfThresh) = 1;
        SignalIdeal = 0;
        for pointer = 2:Headers(X,Y,2)
            BelowThresh = 0;
            for Thresh = 1:NumOfThresh
                if PixelatedBoard(X,Y,4,pointer) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                    BelowThresh = Thresh;
                    break
                end
            end
            if BelowThresh ~= 0
                AboveThresh = 0;
                for Thresh = BelowThresh:NumOfThresh
                    if PixelatedBoard(X,Y,2,pointer) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                        AboveThresh = Thresh;
                    else
                        break
                    end
                end
                if AboveThresh ~= 0
                        NoCSCA_STD_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_STD_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                        %%PCS CODE SECTION
                        if PCS
                            if BelowThresh == 1;
	                        PCS_closetime = PixelatedBoard(X,Y,1,pointer)   + PCS_INTEGRATIONTIME;
                            end
                            if PixelatedBoard(X,Y,1,pointer) <= PCS_closetime;
                                NoCSCA_PCS_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_PCS_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                            end
                        end
                        %%PCS CODE SECTION
                        
                        %%DLR CODE SECTION
                        if DLR
                                if PixelatedBoard(X,Y,1,pointer) <= DLR_closetime;
                                    DLR_HighestThreshReached = max(AboveThresh,DLR_HighestThreshReached);
                                elseif PixelatedBoard(X,Y,1,pointer) <= (DLR_closetime + DLR_WriteTime)
                                else
                                    NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) = NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) + 1;
                                    DLR_HighestThreshReached = AboveThresh;
                                    DLR_closetime = PixelatedBoard(X,Y,1,pointer)   + DLR_INTEGRATIONTIME;
                                end
                           
                        end
                        %%DLR CODE SECTION
                        
                        %%SR CODE SECTION
                        if SR
                            if PixelatedBoard(X,Y,1,pointer) <= SR_closetime
                                SR_Buffer(BelowThresh:AboveThresh) = 1;
                            elseif  PixelatedBoard(X,Y,1,pointer) <= (SR_closetime + SR_WriteTime)
                            else
                                NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) = NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) + 1; %Read out the SR Buffer
                                SR_Buffer(:) = 0;
                                SR_Buffer(BelowThresh:AboveThresh) = 1;
                                SR_closetime = PixelatedBoard(X,Y,1,pointer)   + SR_INTEGRATIONTIME;
                            end
                        end
                        %%SR CODE SECTION
               end
            end
                                    %%IDEAL CODE SECTION
                        if IDEAL
                            if PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1) < 10^-9
                                SignalIdeal = SignalIdeal + PixelatedBoard(X,Y,2,pointer);
                            else
                                AboveThresh = 0;
                                for Thresh = 1:NumOfThresh
                                    if SignalIdeal > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                        AboveThresh = Thresh;
                                    else
                                        break
                                    end
                                end
                                if AboveThresh ~= 0
                                    NoCSCA_IDEAL_Counters(X,Y,1:AboveThresh) = NoCSCA_IDEAL_Counters(X,Y,1:AboveThresh) + 1;
                                end
                                SignalIdeal = 0;
                            end
                        end
                        %%IDEAL CODE SECTION
            NoCSCA_EventsProcessed = NoCSCA_EventsProcessed + 1;
            if rem(NoCSCA_EventsProcessed,Sz/1000) == 0
                NoCSCA_PercentageCompleted = NoCSCA_EventsProcessed/(Sz/100)
            end
        end
        if DLR
            NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) = NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) + 1; %This line ensures the DLR buffer is emptied at the end of the simulation, preventing dropped events
        end

        if SR
            NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) = NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) + 1; %This line ensures the SR buffer is emptied at the end of the simulation, preventing dropped events
        end
    end
end

end




for NoiseOff_CodeFor3x3DyAndHybrid = 1:1
%% Apply CSCAs - 3x3 and Hybrid (Dynamic)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add3x3Dy || Sub3x3Dy || AddHybridDy || SubHybridDy
for TestWhichCSCAs = 1:1
    if Add3x3Dy
        Add3x3Dy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
    if Sub3x3Dy
        Sub3x3Dy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
    if AddHybridDy
        AddHybridDy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
    if SubHybridDy
        SubHybridDy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.
for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixeles to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);
  
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)); %This is the decay corrected signal background on which this event sits

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
                SigB4Array_3x3(1:3,1:3) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_3x3Neighbourhood = zeros(3,3);
                for dx = -1:1
                    for dy = -1:1
                        tempX = X + dx;
                        tempY = Y + dy;
                        if tempX >= 1 && tempX <= XMax && tempY >= 1 && tempY <= YMax
                            localX = 2+dx;
                            localY = 2+dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);
                            SigB4Array_3x3(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent);
                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_3x3Neighbourhood(localX,localY) = max(Signals_3x3Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn));
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
Headers(tempX,tempY,1) = Headers(tempX,tempY,1) + dn;
                            end
                        end
                    end
                end


%%IF YOU WANT TO ADD NOISE, THIS WOULD BE A GOOD POINT TO INPUT IT HERE FOR
%%THE 3X3 AND HYBRID CSCAS

                if Add3x3Dy %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_3x3Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_3x3Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_3x3Neighbourhood));
                    SigB4 = SigB4Array_3x3(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - 2;
                            AY = Y + ddy -2;
                            Add3x3Dy_Counters(AX,AY,BelowThresh:AboveThresh) = Add3x3Dy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub3x3Dy
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:3
                        for probey = 1:3
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_3x3(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_3x3(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_3x3Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == 2 && probey == 2
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub3x3Dy_Counters(X,Y,:) = Sub3x3Dy_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            


            if AddHybridDy
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the maximum neighbourhood
                [SubPixSum, MaxSumIndex] = max(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MaxSumIndex);
                
                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X -2 + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y -2 + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        AddHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) = AddHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end


            if SubHybridDy
                %First determine the smallest signals within the 4 sub neighbourhoods
                SubPixSum = 1000000000; %This line lmits to TeV. This should be enough in almost all instances as an upper limit, but be aware for some nuclear applications perhaps.

                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = SubPixSum;
                
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the minimum neighbourhood score
                [SubPixSum, MinSumIndex] = min(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MinSumIndex);

                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = -SubPixSum;
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X -2 + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y -2 + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        SubHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) = SubHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end
        else
            if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end
        end
    end

    Hybrid_AndOr3x3Dy_EventsProcessed = Hybrid_AndOr3x3Dy_EventsProcessed + 1;
    if rem(Hybrid_AndOr3x3Dy_EventsProcessed,Sz/1000) == 0
        Hybrid_AndOr3x3Dy_PercentageCompleted = Hybrid_AndOr3x3Dy_EventsProcessed/(Sz/100)
    end
end
end

end

for NoiseOff_CodeFor2x2Dy= 1:1
%%Apply CSCAs - 2x2 (Dynamic)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add2x2Dy || Sub2x2Dy
    for TestWhichCSCAs = 1:1
    if Add2x2Dy
        Add2x2Dy_Counters = zeros(XMax, YMax,NumOfThresh);
        Both2x2Dy_EventsProcessed = 0;
    end
    if Sub2x2Dy
        Sub2x2Dy_Counters =zeros(XMax, YMax,NumOfThresh);
        Both2x2Dy_EventsProcessed = 0;
    end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.

    for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixels to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)); %This is the decay corrected signal background on which this event sits

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
            SigB4Array_2x2(1:2,1:2) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_2x2Neighbourhood = zeros(2,2);
                
                for dx = 0:1
                    for dy = 0:1
                        tempX = X+dx;
                        tempY = Y+dy;
                        if tempX <= XMax && tempY <= YMax
                            localX = 1 + dx;
                            localY = 1 + dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);
                            SigB4Array_2x2(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent);
                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_2x2Neighbourhood(localX,localY) = max(Signals_2x2Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn)); %max not sum as you are adding signals which are already decay corrected so include each other... essentially this is a peak track and hold.
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)+dn),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
                                Headers(tempX,tempY,1) = Headers(tempX,tempY,1) + dn;
                            end
                        end
                    end
                end

                if Add2x2Dy %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_2x2Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_2x2Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_2x2Neighbourhood));
                    SigB4 = SigB4Array_2x2(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - 1;
                            AY = Y + ddy -1;
                            Add2x2Dy_Counters(AX,AY,BelowThresh:AboveThresh) = Add2x2Dy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub2x2Dy
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:2
                        for probey = 1:2
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_2x2(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_2x2(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_2x2Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == 1 && probey == 1
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub2x2Dy_Counters(X,Y,:) = Sub2x2Dy_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            else
            if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end
        end
    end

    Both2x2Dy_EventsProcessed = Both2x2Dy_EventsProcessed+1;
    if rem(Both2x2Dy_EventsProcessed,Sz/1000) == 0
        Both2x2Dy_PercentageCompleted = Both2x2Dy_EventsProcessed/(Sz/100)
    end
end
end

end

for NoiseOff_CodeFor3x3StAndHybrid = 1:1
%% Apply CSCAs - 3x3 and Hybrid (Static)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add3x3St || Sub3x3St || AddHybridSt || SubHybridSt
for TestWhichCSCAs = 1:1
    if Add3x3St
        Add3x3St_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
    if Sub3x3St
        Sub3x3St_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
    if AddHybridSt
        AddHybridSt_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
    if SubHybridSt
        SubHybridSt_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.
for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixeles to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)); %This is the decay corrected signal background on which this event sits

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
                SigB4Array_3x3(1:3,1:3) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_3x3Neighbourhood = zeros(3,3);

                switch rem(X-1,3)
                    case 0
                        lx = 0;
                        ux = 2;
                        locx = 1;
                    case 1
                        lx = -1;
                        ux = 1;
                        locx = 2;
                    case 2
                        lx = -2;
                        ux = 0;
                        locx = 3;
                end

                switch rem(Y-1,3)
                    case 0
                        ly = 0;
                        uy = 2;
                        locy = 1;
                    case 1
                        ly = -1;
                        uy = 1;
                        locy = 2;
                    case 2
                        ly = -2;
                        uy = 0;
                        locy = 3;
                end
                
                for dx = lx:ux
                    for dy = ly:uy
                        tempX = X + dx;
                        tempY = Y + dy;
                        if tempX >= 1 && tempX <= XMax && tempY >= 1 && tempY <= YMax
                            localX = locx+dx;
                            localY = locy+dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);
                            SigB4Array_3x3(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent);
                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_3x3Neighbourhood(localX,localY) = max(Signals_3x3Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn));
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
Headers(tempX,tempY,1) = Headers(tempX,tempY,1) + dn;
                            end
                        end
                    end
                end

                if Add3x3St %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_3x3Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_3x3Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_3x3Neighbourhood));
                    SigB4 = SigB4Array_3x3(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - locx;
                            AY = Y + ddy - locy;
                            Add3x3St_Counters(AX,AY,BelowThresh:AboveThresh) = Add3x3St_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub3x3St
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:3
                        for probey = 1:3
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_3x3(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_3x3(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_3x3Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == locx && probey == locy
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub3x3St_Counters(X,Y,:) = Sub3x3St_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            



            if AddHybridSt
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the maximum neighbourhood
                [SubPixSum, MaxSumIndex] = max(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MaxSumIndex);
                
                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X - locx + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y - locy + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        AddHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) = AddHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end


            if SubHybridSt
                %First determine the smallest signals within the 4 sub neighbourhoods
                SubPixSum = 1000000000; %This line lmits to TeV. This should be enough in almost all instances as an upper limit, but be aware for some nuclear applications perhaps.
                
                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = SubPixSum;
                
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the minimum neighbourhood score
                [SubPixSum, MinSumIndex] = min(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MinSumIndex);

                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = - SubPixSum;
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X - locx + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y - locy + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        SubHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) = SubHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end
        else
            if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end
        end
    end

    Hybrid_AndOr3x3St_EventsProcessed = Hybrid_AndOr3x3St_EventsProcessed+1;
    if rem(Hybrid_AndOr3x3St_EventsProcessed,Sz/1000) == 0
        Hybrid_AndOr3x3St_PercentageCompleted = Hybrid_AndOr3x3St_EventsProcessed/(Sz/100)
    end
end
end

end

for NoiseOff_CodeFor2x2St = 1:1
%%Apply CSCAs - 2x2 (Static)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add2x2St || Sub2x2St
    for TestWhichCSCAs = 1:1
    if Add2x2St
        Add2x2St_Counters = zeros(XMax, YMax,NumOfThresh);
        Both2x2St_EventsProcessed = 0;
    end
    if Sub2x2St
        Sub2x2St_Counters = zeros(XMax, YMax,NumOfThresh);
   Both2x2Dy_EventsProcessed = 0;
   end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.
for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixeles to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)); %This is the decay corrected signal background on which this event sits

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
                SigB4Array_2x2(1:2,1:2) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_2x2Neighbourhood= zeros(2,2);
                
                switch rem(X-1,2)
                    case 0
                    lx = 0;
                    ux = 1;
                    locx = 1;
                    case 1
                    lx = -1;
                    ux = 0;
                    locx = 2;
                end

                switch rem(Y-1,2)
                    case 0
                    ly = 0;
                    uy = 1;
                    locy = 1;
                    case 1
                    ly = -1;
                    uy = 0;
                    locy = 2;
                end

                %%%%%SearchingWindow Stuff added
                for dx = lx:ux
                    for dy = ly:uy
                        tempX = X+dx;
                        tempY = Y+dy;
                        if tempX <= XMax && tempY <= YMax
                            localX = locx + dx;
                            localY = locy +dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);
                            SigB4Array_2x2(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent);
                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_2x2Neighbourhood(localX,localY) = max(Signals_2x2Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn));
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
                                Headers(tempX,tempY,1) = Headers(tempX,tempY,1) +dn; %%%%%%%%%%%%%MOVED%%%%%%%%%%%
                            end
                        end
                    end
                end

                if Add2x2St %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_2x2Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_2x2Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_2x2Neighbourhood));
                    SigB4 = SigB4Array_2x2(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - locx;
                            AY = Y + ddy -locy;
                            Add2x2St_Counters(AX,AY,BelowThresh:AboveThresh) = Add2x2St_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub2x2St
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:2
                        for probey = 1:2
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_2x2(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_2x2(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_2x2Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == locx && probey == locy
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub2x2St_Counters(X,Y,:) = Sub2x2St_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            else %%%%%%%%%%%%%%%%ADDED%%%%%%%%%%%%%%%%%
       if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end %%%%%%%%%%%%%%ADDED%%%%%%%%%%%%
        end
    end

    Both2x2St_EventsProcessed = Both2x2St_EventsProcessed + 1;
    if rem(Both2x2St_EventsProcessed,Sz/1000) == 0
        Both2x2St_PercentageCompleted = Both2x2St_EventsProcessed/(Sz/100)
    end
end
end

end



end
end
for NoiseOn = 1:1
if NOISEON == 1

    
for NoiseOn_ApplyFDR = 1:1    
%% Apply FDR counting scheme if requested - In current form, this needs to be implemented before decay correction BUT ALSO IS NOT AFFECTED BY PULSE SHAPE IN THAT CASE
fdrstart = datetime;
for Applying_FDR_ACTs = 1
    if FDR == 1
        NoCSCA_FDR_Counters = zeros(XMax, YMax,NumOfThresh);
        NoCSCA_FDR_EventsProcessed = 0;
        NoCSCA_PercentageCompleted = 0;
        NoCSCA_FDR_updateInterval = Sz*0.05;
        NoCSCA_FDR_nextUpdate = NoCSCA_FDR_updateInterval;
        for X = 1: XMax
            for Y = 1:YMax

%%New FDR version starts here                
PrevSignal = 0;
pointer = 2;
FrameSituation = 0;
OutsideOfFrame = 0;
while pointer <= Headers(X,Y,2) %Known events for this pixel. While needed to allow continue to loop over events in the out of frame section
  
    %%Percentage complete tracker
    NoCSCA_FDR_EventsProcessed = NoCSCA_FDR_EventsProcessed + 1;
    if NoCSCA_FDR_EventsProcessed >= NoCSCA_FDR_nextUpdate
fprintf('FDR: %.2f%% completed\n', (NoCSCA_FDR_EventsProcessed/Sz)*100);
NoCSCA_FDR_nextUpdate = NoCSCA_FDR_nextUpdate + NoCSCA_FDR_updateInterval;
    end
    %%Percentage complete tracker

    if FrameSituation  == 0
        %DO REGULAR SIG CALC
        if OutsideOfFrame == 1
            deltaT = (PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1));
            PrevSignal = CurrentSignal*exp(-DecayConstant*deltaT); %calculate and record residual signal in pixel before current event
        end
        CurrentSignal = PixelatedBoard(X,Y,2,pointer) + PrevSignal; %add current event onto residual signal
        
        %NOISE: Include Noise
            PrevSignal_Noised = PrevSignal + PixelatedBoard(X,Y,5,pointer-1);
            CurrentSignal_Noised = CurrentSignal + PixelatedBoard(X,Y,5,pointer);;
        %COMP TO THRESHOLDS
        BelowThresh = 0;
        for Thresh = 1:NumOfThresh
            if PrevSignal_Noised <= Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                BelowThresh = Thresh;
                break
            end
        end
        if BelowThresh ~= 0
            AboveThresh = 0;
            for Thresh = BelowThresh:NumOfThresh
                if CurrentSignal_Noised > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                    AboveThresh = Thresh;
                else
                    break
                end
            end
            %IF THRESHOLD TRIGGERED
            if AboveThresh ~= 0
                %ADD TO COUNTERS
                NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                %CALC NEW TIMES
                ResetStartTime = PixelatedBoard(X,Y,1,pointer) + FDR_RESETTIME;
                ResetEndTime = ResetStartTime + FDR_RESETDURATION;
                %CHANGE FRAMESITUATION TO 1
                FrameSituation = 1;
            end
        end
        OutsideOfFrame = 1;
    else %if FrameSituation ~= 0, so In Frame
        %CALCULATE DELTAT
        deltaT = (PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1));
        %COMPARE DELTAT TO TIMES TO DETERMINE SITUATIONS
        %IF IN FRAME
        if PixelatedBoard(X,Y,1,pointer) < ResetStartTime
            PrevSignal = CurrentSignal*exp(-DecayConstant*deltaT); %calculate and record residual signal in pixel before current event
            %KEEP COUNTING AS USUAL
            %DO REGULAR SIG CALC
            CurrentSignal = PixelatedBoard(X,Y,2,pointer) + PrevSignal; %add current event onto residual signal

            %NOISE: Include Noise
            PrevSignal_Noised = PrevSignal + PixelatedBoard(X,Y,5,pointer-1);
            CurrentSignal_Noised = CurrentSignal + PixelatedBoard(X,Y,5,pointer);
            %COMP TO THRESHOLDS
            BelowThresh = 0;
            for Thresh = 1:NumOfThresh
                if PrevSignal_Noised <= Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                    BelowThresh = Thresh;
                    break
                end
            end
            if BelowThresh ~= 0
                AboveThresh = 0;
                for Thresh = BelowThresh:NumOfThresh
                    if CurrentSignal_Noised > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                        AboveThresh = Thresh;
                    else
                        break
                    end
                end
                %IF THRESHOLD TRIGGERED
                if AboveThresh ~= 0
                    %ADD TO COUNTERS
                    NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_FDR_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                end
            end

            %ELESIF IN RESET
        elseif PixelatedBoard(X,Y,1,pointer) > ResetEndTime %(IMPLICIT BEYOND FRAME, also implicitly ignores events in reset period)
            %PREVSIGNAL = 0
            PrevSignal = 0;
            %FRAMESITUATION == 0
            FrameSituation = 0;
            OutsideOfFrame = 0;
            %BREAK
            continue
            %END
        end
    end
    pointer = pointer +1;
end
%%New FDR version ends
            end
        end
    end
end
fprintf('FDR: 100%% completed\n');

%THIS MODULE SHOULD BE TESTED IN A SANDBOX
  fdrstop = datetime;
  fdrduration = fdrstop - fdrstart

end

for NoiseOn_DecayCorrectSignals = 1:1
%% Decay correct signals to be a running total
switch non_impulsePulseShape
case 1
    
case 2
for X = 1: XMax
    for Y = 1:YMax
        for pointer = 2:Headers(X,Y,2)
            deltaT = PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1); %calculate time since last event
            PixelatedBoard(X,Y,4,pointer) = PixelatedBoard(X,Y,2,pointer-1)*exp(-DecayConstant*deltaT); %calculate and record residual signal in pixel before current event
            PixelatedBoard(X,Y,2,pointer) = PixelatedBoard(X,Y,2,pointer) + PixelatedBoard(X,Y,4,pointer); %add current event onto residual signal
        end
    end
end

if DEBUGGING == 1
fprintf('Saving Output Data For Debug Logs...');
startA = datetime;
save(strcat(LABEL,'_DBG4'), 'PixelatedBoard', '-v7.3')
stopB = datetime;
fprintf('Completed in (hr:min:sec): %s\n', stopB-startA);
end

case 3 
%%%%%%%%%%%%%%%%%%%%%%%%%
%3). From empirical data (slowest, empirical model)\n',[1,2,3],REPEATEDINCORRECTCHOICEERRORMESSAGE);
%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

A = datetime;

for NoiseOn_ApplyNonFDR_ACTS = 1:1
%% Apply ACTS
%Define Thresholds
algstart = datetime;
Thresholds = THRESHOLDS;
NumOfThresh = max(size(Thresholds));
% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
NoCSCA_STD_Counters = zeros(XMax, YMax,NumOfThresh);
if PCS
    NoCSCA_PCS_Counters = zeros(XMax, YMax,NumOfThresh);
    PCS_closetime = 0;
end
if DLR
    NoCSCA_DLR_Counters = -1 * ones(XMax, YMax,NumOfThresh); %-1 as initialising the buffer causes a single count to go to each threshold, so this means the system activation only will return 0 counts
end
if SR
    NoCSCA_SR_Counters = -1 * ones(XMax, YMax,NumOfThresh); %-1 as initialising the buffer causes a single count to go to each threshold, so this means the system activation only will return 0 counts
end
if IDEAL
NoCSCA_IDEAL_Counters = zeros(XMax, YMax,NumOfThresh);
end

NoCSCA_EventsProcessed = 0;
for X = 1:XMax
    for Y = 1:YMax
        DLR_HighestThreshReached = NumOfThresh;
        DLR_closetime = - abs(2.0*PixelatedBoard(X,Y,1,2) + DLR_INTEGRATIONTIME);
        SR_closetime = - abs(2.0*PixelatedBoard(X,Y,1,2) + SR_INTEGRATIONTIME);
        SR_Buffer(1:NumOfThresh) = 1;
        SignalIdeal = 0;
        for pointer = 2:Headers(X,Y,2)
            BelowThresh = 0;
            for Thresh = 1:NumOfThresh
                if (PixelatedBoard(X,Y,4,pointer)+PixelatedBoard(X,Y,5,pointer-1)) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                    BelowThresh = Thresh;
                    break
                end
            end
            if BelowThresh ~= 0
                AboveThresh = 0;
                for Thresh = BelowThresh:NumOfThresh
                    if (PixelatedBoard(X,Y,2,pointer)+PixelatedBoard(X,Y,5,pointer)) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                        AboveThresh = Thresh;
                    else
                        break
                    end
                end
                if AboveThresh ~= 0
                        NoCSCA_STD_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_STD_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                        %%PCS CODE SECTIONNOTE
                        if PCS
                            if BelowThresh == 1;
	                        PCS_closetime = PixelatedBoard(X,Y,1,pointer)   + PCS_INTEGRATIONTIME;
                            end
                            if PixelatedBoard(X,Y,1,pointer) <= PCS_closetime;
                                NoCSCA_PCS_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_PCS_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                            end
                        end
                        %%PCS CODE SECTION
                        
                        %%DLR CODE SECTION
                        if DLR
                                if PixelatedBoard(X,Y,1,pointer) <= DLR_closetime;
                                    DLR_HighestThreshReached = max(AboveThresh,DLR_HighestThreshReached);
                                elseif PixelatedBoard(X,Y,1,pointer) <= (DLR_closetime + DLR_WriteTime)
                                else
                                    NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) = NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) + 1;
                                    DLR_HighestThreshReached = AboveThresh;
                                    DLR_closetime = PixelatedBoard(X,Y,1,pointer)   + DLR_INTEGRATIONTIME;
                                end
                           
                        end
                        %%DLR CODE SECTION
                        
                        %%SR CODE SECTION
                        if SR
                            if PixelatedBoard(X,Y,1,pointer) <= SR_closetime
                                SR_Buffer(BelowThresh:AboveThresh) = 1;
                            elseif  PixelatedBoard(X,Y,1,pointer) <= (SR_closetime + SR_WriteTime)
                            else
                                NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) = NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) + 1; %Read out the SR Buffer
                                SR_Buffer(:) = 0;
                                SR_Buffer(BelowThresh:AboveThresh) = 1;
                                SR_closetime = PixelatedBoard(X,Y,1,pointer)   + SR_INTEGRATIONTIME;
                            end
                        end
                        %%SR CODE SECTION
               end
            end
                                    %%IDEAL CODE SECTION
                        if IDEAL
                            if PixelatedBoard(X,Y,1,pointer) - PixelatedBoard(X,Y,1,pointer-1) < 10^-9
                                SignalIdeal = SignalIdeal + PixelatedBoard(X,Y,2,pointer);
                            else
                                AboveThresh = 0;
                                for Thresh = 1:NumOfThresh
                                    if SignalIdeal > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                        AboveThresh = Thresh;
                                    else
                                        break
                                    end
                                end
                                if AboveThresh ~= 0
                                    NoCSCA_IDEAL_Counters(X,Y,1:AboveThresh) = NoCSCA_IDEAL_Counters(X,Y,1:AboveThresh) + 1;
                                end
                                SignalIdeal = 0;
                            end
                        end
                        %%IDEAL CODE SECTION
            NoCSCA_EventsProcessed = NoCSCA_EventsProcessed + 1;
            if rem(NoCSCA_EventsProcessed,Sz/1000) == 0
                NoCSCA_PercentageCompleted = NoCSCA_EventsProcessed/(Sz/100)
            end
        end
        if DLR
            NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) = NoCSCA_DLR_Counters(X,Y,1:DLR_HighestThreshReached) + 1; %This line ensures the DLR buffer is emptied at the end of the simulation, preventing dropped events
        end

        if SR
            NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) = NoCSCA_SR_Counters(X,Y,1:sum(SR_Buffer)) + 1; %This line ensures the SR buffer is emptied at the end of the simulation, preventing dropped events
        end
    end
end


end



for NoiseOn_CodeFor3x3DyAndHybrid = 1:1    
%% Apply CSCAs - 3x3 and Hybrid (Dynamic)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add3x3Dy || Sub3x3Dy || AddHybridDy || SubHybridDy
for TestWhichCSCAs = 1:1
    if Add3x3Dy
        Add3x3Dy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
    if Sub3x3Dy
        Sub3x3Dy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
    if AddHybridDy
        AddHybridDy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
    if SubHybridDy
        SubHybridDy_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3Dy_EventsProcessed = 0;
    end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.
for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixeles to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);
        
        %%NOISE: modified
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)-1); %This is the decay corrected signal background on which this event sits
        %%NOISE: modified

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end 

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
                SigB4Array_3x3(1:3,1:3) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_3x3Neighbourhood = zeros(3,3);
                for dx = -1:1
                    for dy = -1:1
                        tempX = X + dx;
                        tempY = Y + dy;
                        if tempX >= 1 && tempX <= XMax && tempY >= 1 && tempY <= YMax
                            localX = 2+dx;
                            localY = 2+dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);
                            
                            %%ADD NOISE - local start
                            PrevNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1)-1);
                            %ADD NOISE - local end
                            
                            %%NOISE: modified
                            SigB4Array_3x3(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent) + PrevNoise; %Line modified for noise
                            %%NOISE: modified

                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_3x3Neighbourhood(localX,localY) = max(Signals_3x3Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn));
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end

                                %%ADD NOISE - local start
                                NextNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1));
                                ExtrapolatedNoise = ( PrevNoise*TSinceLastEvent + NextNoise*(NoiseFreq-TSinceLastEvent) ) / sqrt( (TSinceLastEvent^2) + (NoiseFreq-TSinceLastEvent)^2);
                                Signals_3x3Neighbourhood(localX,localY) = Signals_3x3Neighbourhood(localX,localY) + ExtrapolatedNoise;
                                %ADD NOISE - local end

                                Headers(tempX,tempY,1) = Headers(tempX,tempY,1) + dn;
                            end
                        end
                    end
                end


%%IF YOU WANT TO ADD NOISE, THIS WOULD BE A GOOD POINT TO INPUT IT HERE FOR
%%THE 3X3 AND HYBRID CSCAS

                if Add3x3Dy %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_3x3Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_3x3Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_3x3Neighbourhood));
                    SigB4 = SigB4Array_3x3(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - 2;
                            AY = Y + ddy -2;
                            Add3x3Dy_Counters(AX,AY,BelowThresh:AboveThresh) = Add3x3Dy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub3x3Dy
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:3
                        for probey = 1:3
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_3x3(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_3x3(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_3x3Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == 2 && probey == 2
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub3x3Dy_Counters(X,Y,:) = Sub3x3Dy_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            


            if AddHybridDy
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the maximum neighbourhood
                [SubPixSum, MaxSumIndex] = max(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MaxSumIndex);
                
                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X -2 + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y -2 + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        AddHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) = AddHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end


            if SubHybridDy
                %First determine the smallest signals within the 4 sub neighbourhoods
                SubPixSum = 1000000000; %This line lmits to TeV. This should be enough in almost all instances as an upper limit, but be aware for some nuclear applications perhaps.

                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = SubPixSum;
                
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the minimum neighbourhood score
                [SubPixSum, MinSumIndex] = min(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MinSumIndex);

                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = -SubPixSum;
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X -2 + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y -2 + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        SubHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) = SubHybridDy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end
        else
            if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end
        end
    end

    Hybrid_AndOr3x3Dy_EventsProcessed = Hybrid_AndOr3x3Dy_EventsProcessed + 1;
    if rem(Hybrid_AndOr3x3Dy_EventsProcessed,Sz/1000) == 0
        Hybrid_AndOr3x3Dy_PercentageCompleted = Hybrid_AndOr3x3Dy_EventsProcessed/(Sz/100)
    end
end
end

end

for NoiseOn_CodeFor2x2Dy = 1:1    
%%Apply CSCAs - 2x2 (Dynamic)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add2x2Dy || Sub2x2Dy
    for TestWhichCSCAs = 1:1
    if Add2x2Dy
        Add2x2Dy_Counters = zeros(XMax, YMax,NumOfThresh);
        Both2x2Dy_EventsProcessed = 0;
    end
    if Sub2x2Dy
        Sub2x2Dy_Counters =zeros(XMax, YMax,NumOfThresh);
        Both2x2Dy_EventsProcessed = 0;
    end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.

    for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixels to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);

        %%NOISE: modified
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)-1); %This is the decay corrected signal background on which this event sits
        %%NOISE: modified

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
            SigB4Array_2x2(1:2,1:2) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_2x2Neighbourhood = zeros(2,2);
                
                for dx = 0:1
                    for dy = 0:1
                        tempX = X+dx;
                        tempY = Y+dy;
                        if tempX <= XMax && tempY <= YMax
                            localX = 1 + dx;
                            localY = 1 + dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);

                            %%ADD NOISE - local start
                            PrevNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1)-1);
                            %ADD NOISE - local end

                            %%NOISE: modified
                            SigB4Array_2x2(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent) + PrevNoise; %Line modified for noise
                            %%NOISE: modified

                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_2x2Neighbourhood(localX,localY) = max(Signals_2x2Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn)); %max not sum as you are adding signals which are already decay corrected so include each other... essentially this is a peak track and hold.
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)+dn),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
                                
                                %%ADD NOISE - local start
                                NextNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1));
                                ExtrapolatedNoise = ( PrevNoise*TSinceLastEvent + NextNoise*(NoiseFreq-TSinceLastEvent) ) / sqrt( (TSinceLastEvent^2) + (NoiseFreq-TSinceLastEvent)^2);
                                Signals_2x2Neighbourhood(localX,localY) = Signals_2x2Neighbourhood(localX,localY) + ExtrapolatedNoise;
                                %ADD NOISE - local end

                                Headers(tempX,tempY,1) = Headers(tempX,tempY,1) + dn;
                            end
                        end
                    end
                end

                if Add2x2Dy %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_2x2Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_2x2Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_2x2Neighbourhood));
                    SigB4 = SigB4Array_2x2(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - 1;
                            AY = Y + ddy -1;
                            Add2x2Dy_Counters(AX,AY,BelowThresh:AboveThresh) = Add2x2Dy_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub2x2Dy
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:2
                        for probey = 1:2
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_2x2(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_2x2(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_2x2Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == 1 && probey == 1
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub2x2Dy_Counters(X,Y,:) = Sub2x2Dy_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            else
            if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end
        end
    end

    Both2x2Dy_EventsProcessed = Both2x2Dy_EventsProcessed+1;
    if rem(Both2x2Dy_EventsProcessed,Sz/1000) == 0
        Both2x2Dy_PercentageCompleted = Both2x2Dy_EventsProcessed/(Sz/100)
    end
end
end


end

for NoiseOn_CodeFor3x3StAndHybrid = 1:1    
%% Apply CSCAs - 3x3 and Hybrid (Static)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add3x3St || Sub3x3St || AddHybridSt || SubHybridSt
for TestWhichCSCAs = 1:1
    if Add3x3St
        Add3x3St_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
    if Sub3x3St
        Sub3x3St_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
    if AddHybridSt
        AddHybridSt_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
    if SubHybridSt
        SubHybridSt_Counters = zeros(XMax, YMax,NumOfThresh);
        Hybrid_AndOr3x3St_EventsProcessed = 0;
    end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.
for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixeles to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);

        %%NOISE: modified
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)-1); %This is the decay corrected signal background on which this event sits
        %%NOISE: modified

        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
                SigB4Array_3x3(1:3,1:3) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_3x3Neighbourhood = zeros(3,3);

                switch rem(X-1,3)
                    case 0
                        lx = 0;
                        ux = 2;
                        locx = 1;
                    case 1
                        lx = -1;
                        ux = 1;
                        locx = 2;
                    case 2
                        lx = -2;
                        ux = 0;
                        locx = 3;
                end

                switch rem(Y-1,3)
                    case 0
                        ly = 0;
                        uy = 2;
                        locy = 1;
                    case 1
                        ly = -1;
                        uy = 1;
                        locy = 2;
                    case 2
                        ly = -2;
                        uy = 0;
                        locy = 3;
                end
                
                for dx = lx:ux
                    for dy = ly:uy
                        tempX = X + dx;
                        tempY = Y + dy;
                        if tempX >= 1 && tempX <= XMax && tempY >= 1 && tempY <= YMax
                            localX = locx+dx;
                            localY = locy+dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);

                            %%ADD NOISE - local start
                            PrevNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1)-1);
                            %ADD NOISE - local end

                            %%NOISE: modified
                            SigB4Array_3x3(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent) + PrevNoise; %Line modified for noise
                            %%NOISE: modified

                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_3x3Neighbourhood(localX,localY) = max(Signals_3x3Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn));
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
                                
                                %%ADD NOISE - local start
                                NextNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1));
                                ExtrapolatedNoise = ( PrevNoise*TSinceLastEvent + NextNoise*(NoiseFreq-TSinceLastEvent) ) / sqrt( (TSinceLastEvent^2) + (NoiseFreq-TSinceLastEvent)^2);
                                Signals_3x3Neighbourhood(localX,localY) = Signals_3x3Neighbourhood(localX,localY) + ExtrapolatedNoise;
                                %ADD NOISE - local end

Headers(tempX,tempY,1) = Headers(tempX,tempY,1) + dn;
                            end
                        end
                    end
                end

                if Add3x3St %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_3x3Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_3x3Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_3x3Neighbourhood));
                    SigB4 = SigB4Array_3x3(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - locx;
                            AY = Y + ddy - locy;
                            Add3x3St_Counters(AX,AY,BelowThresh:AboveThresh) = Add3x3St_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub3x3St
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:3
                        for probey = 1:3
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_3x3(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_3x3(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_3x3Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == locx && probey == locy
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub3x3St_Counters(X,Y,:) = Sub3x3St_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            



            if AddHybridSt
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the maximum neighbourhood
                [SubPixSum, MaxSumIndex] = max(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MaxSumIndex);
                
                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X - locx + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y - locy + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        AddHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) = AddHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end


            if SubHybridSt
                %First determine the smallest signals within the 4 sub neighbourhoods
                SubPixSum = 1000000000; %This line lmits to TeV. This should be enough in almost all instances as an upper limit, but be aware for some nuclear applications perhaps.
                
                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = SubPixSum;
                
                %First determine the sub signals within the 4 sub neighbourhoods
                local_sums = conv2(Signals_3x3Neighbourhood, ones(2), 'valid');
                %Then find the minimum neighbourhood score
                [SubPixSum, MinSumIndex] = min(local_sums(:));
                %Then find the coordinated for the largest bottom left
                %pixel in this array
                [Subx, Suby] = ind2sub(size(local_sums), MinSumIndex);

                %Next find the coordinates for the largest scoring pixel
                %within the determined subneighbourhood
                Signals_3x3Neighbourhood(SigB4Array_3x3 == -317) = - SubPixSum;
                IndPixSig = max(max(Signals_3x3Neighbourhood(Subx:Subx+1, Suby:Suby+1)));
                [ddx, ddy] = find(Signals_3x3Neighbourhood == IndPixSig);
                %Finally, convery these coordinates to those of the global
                %X,Y system.
                
                SigB4 = SigB4Array_3x3(ddx,ddy);
                BelowThresh = 0;
                for Thresh = 1:NumOfThresh
                    if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                        BelowThresh = Thresh;
                        break
                    end
                end
                if BelowThresh ~= 0
                    AboveThresh = 0;
                    for Thresh = BelowThresh:NumOfThresh
                        if SubPixSum > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                            AboveThresh = Thresh;
                        else
                            break
                        end
                    end
                    if AboveThresh ~= 0
                        AX = X - locx + ddx(1); % -2 is to correct from 3x3 coordinated to X (i.e., as it is 2x2 centred, newX = 2 implies originalX + 0)
                        AY = Y - locy + ddy(1); % -2 is to correct from 3x3 coordinated to Y (i.e., as it is 2x2 centred, newY = 2 implies originalY + 0)
                        SubHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) = SubHybridSt_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                    end
                end
            end
        else
            if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end
        end
    end

    Hybrid_AndOr3x3St_EventsProcessed = Hybrid_AndOr3x3St_EventsProcessed+1;
    if rem(Hybrid_AndOr3x3St_EventsProcessed,Sz/1000) == 0
        Hybrid_AndOr3x3St_PercentageCompleted = Hybrid_AndOr3x3St_EventsProcessed/(Sz/100)
    end
end
end



end

for NoiseOn_CodeFor2x2St = 1:1
%%Apply CSCAs - 2x2 (Static)
%%NOTE THAT IF NOISE IS TO BE ADDED TO THE SYSTEM, IT SHOULD BE MODELED
%%INTO EACH PIXEL IN THE CSCA NEIGHBOURHOODS. THE HYBRID ONES NEED TO ADD
%%THE NOISE BEFORE THE NEIGHBOURHOODING.

%%As different CSCA search areas will involve different events being
%%counted, a single pass through and master list will not be viable.
%%Instead, the masterlist will need to be used, and reset, for each of the
%%3x3 and 2x2 cases. A single pass through with Hybrid should still be
%%viable though. For the same reason, Dynamic and Static CSCA will each
%%require their own pass throughs. 

% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
if Add2x2St || Sub2x2St
    for TestWhichCSCAs = 1:1
    if Add2x2St
        Add2x2St_Counters = zeros(XMax, YMax,NumOfThresh);
        Both2x2St_EventsProcessed = 0;
    end
    if Sub2x2St
        Sub2x2St_Counters = zeros(XMax, YMax,NumOfThresh);
   Both2x2Dy_EventsProcessed = 0;
   end
end

%MasterList = MasterList_Backup;
MasterList(:,4) = ones(1,size(MasterList,1));
Headers(:,:,1) = 2;

%These are all currently designed to use the decay corrected signals. For
%impulses, the old CSCA codes could be used instead.
for nn = 1:Sz           %Need to step through all events within master list to allow events within adjacent pixeles to be processed in correct order
    if MasterList(nn,4) == 1 %prevents events being repeat counted THIS COULD BE AN ISSUE!
        CSCA_Trigger = 0; %%????????????
        X = MasterList(nn,2);
        Y = MasterList(nn,3);

        %%NOISE: modified
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)); %This is the signal plus decay background within the pixel at this time
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1)) + PixelatedBoard(X,Y,5,Headers(X,Y,1)-1); %This is the decay corrected signal background on which this event sits
        %%NOISE: modified
        
        for Thresh = 1:NumOfThresh %This needs to be done first as ANY threshold triggered locally will result in the CSCA being triggered, but what this means differs between CSCAs and the signal may not yet have ben fully counted, so the above threshold part cannot be assessed.
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
            SearchingWindow = PixelatedBoard(X,Y,1,Headers(X,Y,1)) + CSCA_SearchingTime;
                SigB4Array_2x2(1:2,1:2) = -317;
                %Calculate signals and thresholds breached in each pixel of
                %neighbourhood (fastest when Add and Sub both used)
                Signals_2x2Neighbourhood= zeros(2,2);
                
                switch rem(X-1,2)
                    case 0
                    lx = 0;
                    ux = 1;
                    locx = 1;
                    case 1
                    lx = -1;
                    ux = 0;
                    locx = 2;
                end

                switch rem(Y-1,2)
                    case 0
                    ly = 0;
                    uy = 1;
                    locy = 1;
                    case 1
                    ly = -1;
                    uy = 0;
                    locy = 2;
                end

                %%%%%SearchingWindow Stuff added
                for dx = lx:ux
                    for dy = ly:uy
                        tempX = X+dx;
                        tempY = Y+dy;
                        if tempX <= XMax && tempY <= YMax
                            localX = locx + dx;
                            localY = locy +dy;
                            TSinceLastEvent = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)-1);
                            
                            %%ADD NOISE - local start
                            PrevNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1)-1);
                            %ADD NOISE - local end

                            %%NOISE: modified
                            SigB4Array_2x2(localX,localY) = PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)-1) * exp(-DecayConstant*TSinceLastEvent) + PrevNoise;
                            %%NOISE: modified
                            
                            dn = 0;
                            if PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                while PixelatedBoard(tempX,tempY,1,Headers(tempX,tempY,1)+dn) <= SearchingWindow
                                    Signals_2x2Neighbourhood(localX,localY) = max(Signals_2x2Neighbourhood(localX,localY), PixelatedBoard(tempX,tempY,2,Headers(tempX,tempY,1)+dn));
                                    MasterList(PixelatedBoard(tempX,tempY,3,Headers(tempX,tempY,1)),4) = -31745; %Set unprocessed flag to 0 for this event
                                    dn = dn + 1;
                                end
                                
                                %%ADD NOISE - local start
                                NextNoise = PixelatedBoard(tempX,tempY,5,Headers(tempX,tempY,1));
                                ExtrapolatedNoise = ( PrevNoise*TSinceLastEvent + NextNoise*(NoiseFreq-TSinceLastEvent) ) / sqrt( (TSinceLastEvent^2) + (NoiseFreq-TSinceLastEvent)^2);
                                Signals_2x2Neighbourhood(localX,localY) = Signals_2x2Neighbourhood(localX,localY) + ExtrapolatedNoise;
                                %ADD NOISE - local end

                                Headers(tempX,tempY,1) = Headers(tempX,tempY,1) +dn; %%%%%%%%%%%%%MOVED%%%%%%%%%%%
                            end
                        end
                    end
                end

                if Add2x2St %Note that this implementation of an additive 3x3 is based on the selection the output pixel and silencing of others before threshold, implying that the output signal may not count or may not trigger all relevant thresholds, based solely on what thresholds are available in the assigned write out pixel.
                    % Find the maximum value and its linear index and then
                    % convert this to the coordinates for the X and Y
                    [~, max_idpix] = max(Signals_2x2Neighbourhood(:));
                    [ddx, ddy] = ind2sub(size(Signals_2x2Neighbourhood), max_idpix);
                    SumSig = sum(sum(Signals_2x2Neighbourhood));
                    SigB4 = SigB4Array_2x2(ddx,ddy);
                    BelowThresh = 0;
                    for Thresh = 1:NumOfThresh
                        if SigB4 < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                            BelowThresh = Thresh;
                            break
                        end
                    end
                    if BelowThresh ~= 0
                        AboveThresh = 0;
                        for Thresh = BelowThresh:NumOfThresh
                            if SumSig > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                AboveThresh = Thresh;
                            else
                                break
                            end
                        end
                        if AboveThresh ~= 0
                            AX = X + ddx - locx;
                            AY = Y + ddy -locy;
                            Add2x2St_Counters(AX,AY,BelowThresh:AboveThresh) = Add2x2St_Counters(AX,AY,BelowThresh:AboveThresh) + 1;
                        end
                    end
                end

                if Sub2x2St
                    ThresholdsPotentiallyTriggered = zeros(1,1,NumOfThresh);
                    Suppress = 0;
                    for probex = 1:2
                        for probey = 1:2
                            if Suppress < 2
                                BelowThresh = 0;
                                if SigB4Array_2x2(probex,probey) ~= -317 %If pixel exists
                                    for Thresh = 1:NumOfThresh
                                        if SigB4Array_2x2(probex,probey) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                                            BelowThresh = Thresh;
                                            break
                                        end
                                    end
                                    if BelowThresh ~= 0
                                        AboveThresh = 0;
                                        for Thresh = BelowThresh:NumOfThresh
                                            if Signals_2x2Neighbourhood(probex,probey) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                                                AboveThresh = Thresh;
                                            else
                                                break
                                            end
                                        end
                                        if AboveThresh ~= 0
                                            Suppress = Suppress + 1;
                                            if probex == locx && probey == locy
                                                ThresholdsPotentiallyTriggered(1,1,BelowThresh:AboveThresh) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if Suppress < 2
                        Sub2x2St_Counters(X,Y,:) = Sub2x2St_Counters(X,Y,:) + ThresholdsPotentiallyTriggered(1,1,:);
                    end
                end
            else %%%%%%%%%%%%%%%%ADDED%%%%%%%%%%%%%%%%%
       if Headers(X,Y,1) < Headers(X,Y,2)
            Headers(X,Y,1) = Headers(X,Y,1) + 1;
            else
                Headers(X,Y,1) = Headers(X,Y,2);
            end %%%%%%%%%%%%%%ADDED%%%%%%%%%%%%
        end
    end

    Both2x2St_EventsProcessed = Both2x2St_EventsProcessed + 1;
    if rem(Both2x2St_EventsProcessed,Sz/1000) == 0
        Both2x2St_PercentageCompleted = Both2x2St_EventsProcessed/(Sz/100)
    end
end
end


end


end







end

    for DebugCheckAndSave = 1:1
        if DEBUGGING == 1
            varsAfter = who;
            newVars = setdiff(varsAfter, varsBefore);
            outputFileName = strcat(Label,'_DDD_PixBrdAllProccessed');
            save(outputFileName, newVars{:})

        elseif DEBUGGING == 2
            varsAfter = who;
            newVars = setdiff(varsAfter, varsBefore);
            outputFileName = strcat(Label,'_DDD_PixBrdAllProccessed_v7pt3');
            save(outputFileName, newVars{:}, '-v7.3')
        end
    end

end



if SectionsToRun(3) == 1
%FINALLY ADJUST FOR HOT PIXEL NOISE/ ALL PIXEL NOISE INDUCED COUNTS (SEE A4 NOTEBOOK)
ArrayOfCountsToAdd = zeros(XMax,YMax,length(THRESHOLDS));

%NoiseFreq2 value can be adjusted later if needed so that noise can have 2 different frequencies?

for X = 1:XMax
    for Y = 1:YMax

%Determine minimum number of time boxes needed to cover events (equal to
%number of noise samples made)
EventTimesList = PixelatedBoard(X,Y,1,2:end-1);
Sigma = NoiseLevelMap_Adjusted(X,Y);
ArrayOfCountsToAdd(X,Y,:) = CalculateNoiseBasedCounts(EventTimesList, Sigma, SigmaTolerance, THRESHOLDS, NoiseFreq2,t2-t1);

    end
end




end
%%




    for DebugCheckAndSave = 1:1
        if DEBUGGING == 1
            varsAfter = who;
            newVars = setdiff(varsAfter, varsBefore);
            outputFileName = strcat(Label,'_DDD_HotPixCntsCalcuted');
            save(outputFileName, newVars{:})

        elseif DEBUGGING == 2
            varsAfter = who;
            newVars = setdiff(varsAfter, varsBefore);
            outputFileName = strcat(Label,'_DDD_HotPixCntsCalcuted_v7pt3');
            save(outputFileName, newVars{:}, '-v7.3')
        else
            outputFileName = strcat(Label,'_WholeSimFinished_1S',num2str(SectionsToRun(1)),'_2S',num2str(SectionsToRun(2)),'_3S',num2str(SectionsToRun(3)),'v7pt3');
            save(outputFileName, '-v7.3')
        end
    end


end


%AND THEN COPY AND PASTE THE ACTS AND CSCAS?
%THEN FIX NOISE? (CALIBRATION MAP AND PERHAPS WE SHOULDN'T CALIBATE IN MAP AS NOISE IS ALREADY CALIBRATED BY HAVING MEAN OF 0. To CHANGE THIS, ADD NOISE AFTER BASED ON SOME KIND OF 3 SIGMA CALCULATION)?




%ADD SOME NOISE LATER











%{
GenerateCIEMap(ComsolFileAddress)
%


InterptedCIEValues = interpn(xVals, yVals, zVals, MaxSignalStructured3DMap, xqs, yqs, zqs, 'linear', 'extrap');


NEED To ADJUST XVAL, YVAL AND ZVAL TO GO FROM 0 TO MAX DIMENSIONS, NOT
-MAX/2 TO MAX/2
 % Define grid axes based on Vals dimensions
    [xlength, ylength, zlength] = size(vals);
    xGrid = 1:xlength;
    yGrid = 1:ylength;
    zGrid = 1:zlength;

Events = [1 1 1.5];
Events1 = vals(1, 1, 1)
Events2 = vals(1, 1, 2)
EventsInterp = (Events1+Events2)/2

    % Extract query coordinates
    xq = Events(:,1);
    yq = Events(:,2);
    zq = Events(:,3);

    % Perform trilinear interpolation
    interpolatedVals = interpn(xGrid, yGrid, zGrid, vals, xq, yq, zq, 'linear')





    % Extract query coordinates AND NEED TO CONVERT TO NUMERICAL FORM OF
    % MAP. E.G. normalisedx/xstepsize+1
    xq = Events(:,1);
    yq = Events(:,2);
    zq = Events(:,3);

    % Perform trilinear interpolation
    interpolatedVals = interpn(xVals, yVals, zVals, Structured4DMap_AllTimes, xq, yq, zq, 'linear')



%}



