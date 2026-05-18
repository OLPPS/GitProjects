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

function CoGI_v3_1_00_CERN_RunOne(ParamsFile,GeometryFile,NoiseParametersFile,MCDataMatFileName,t1,t2, RANDOMSEED,CIEtolerance, SectionsToRun, MapAddresses, SigmaTolerance, VERBOSE,DEBUGGING,Label,THRESHOLDS)

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
    %clear GateDataFile xpixIDs ypixIDs ValidPixels


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

%%
ss = 0.25
    bins = [0:ss:150,inf];
    RawCodes = EventsToTransfer & (CScase_plusCS == 5)';
    CSEventCodes = EventsToTransfer & (CScase_plusCS ~= 5)';

    CntsOfArrived = SpecAnalysis(GateDataFile_plusCS(RawCodes,energy).*1000,1,bins);
    CntsOfPlusCS_CSonly = SpecAnalysis(GateDataFile_plusCS(CSEventCodes,energy).*CIE_Values(CSEventCodes).*1000,1,bins);
    CntsOfCIE = SpecAnalysis(GateDataFile(RawCodes,energy).*CIE_Values(RawCodes).*1000,1,bins);
    CntsOfPlusCS_CIE = SpecAnalysis(GateDataFile_plusCS(EventsToTransfer,energy).*CIE_Values(EventsToTransfer).*1000,1,bins);

    archivedCntsOfPileup = zeros(1,length(bins));
    archivedCntsOfPileup_PlusCS = zeros(1,length(bins));
    archivedCntsOfPileup_PlusCS_plusCIE = zeros(1,length(bins));

for X = 1:PixelsInXDirection
    for Y = 1: PixelsInYDirection
    PileupPix = RawCodes & GateDataFile_plusCS(:,xpixnum) == X & GateDataFile_plusCS(:,ypixnum) == Y;
    ww = 100*10^(-9);
    RoughBinnedData = roughBin(GateDataFile_plusCS(PileupPix,energy).*1000,GateDataFile_plusCS(PileupPix,time),ww);
    CntsOfPileup = SpecAnalysis(RoughBinnedData,1,bins);
    archivedCntsOfPileup(:) = archivedCntsOfPileup(:) + CntsOfPileup(:);
    
    PileupPix = EventsToTransfer & GateDataFile_plusCS(:,xpixnum) == X & GateDataFile_plusCS(:,ypixnum) == Y;
    RoughBinnedData_plusCS = roughBin(GateDataFile_plusCS(PileupPix,energy).*1000,GateDataFile_plusCS(PileupPix,time),ww);
    CntsOfPileup_PlusCS = SpecAnalysis(RoughBinnedData_plusCS,1,bins);
    archivedCntsOfPileup_PlusCS = archivedCntsOfPileup_PlusCS + CntsOfPileup_PlusCS;
    
    RoughBinnedData_plusCS_plusCIE = roughBin(GateDataFile_plusCS(PileupPix,energy).*CIE_Values(PileupPix).*1000,GateDataFile_plusCS(PileupPix,time),ww);
    CntsOfPileup_PlusCS_plusCIE = SpecAnalysis(RoughBinnedData_plusCS_plusCIE,1,bins);
    archivedCntsOfPileup_PlusCS_plusCIE = archivedCntsOfPileup_PlusCS_plusCIE + CntsOfPileup_PlusCS_plusCIE;
    



    end
end
    
CntsOfPileup = archivedCntsOfPileup;
CntsOfPileup_PlusCS = archivedCntsOfPileup_PlusCS;
CntsOfPileup_PlusCS_PlusCIE = archivedCntsOfPileup_PlusCS_plusCIE;
    

    RoughBinnedData2 = roughBin(GateDataFile_plusCS(RawCodes,energy).*1000.*CIE_Values(RawCodes),GateDataFile_plusCS(RawCodes,time),ww);
    CntsOfPileup_CIE = SpecAnalysis(RoughBinnedData2,1,bins);
    RoughBinnedData_PlusCS2 = roughBin(GateDataFile_plusCS(EventsToTransfer,energy).*CIE_Values(EventsToTransfer).*1000,GateDataFile_plusCS(EventsToTransfer,time),ww);
    CntsOfPileup_PlusCS_PluseCIE = SpecAnalysis(RoughBinnedData_PlusCS2,1,bins);




    figure(1)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfArrived')

    figure(2)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfPlusCS')

    figure(3)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfCIE')

    figure(4)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfPileup')

    figure(5)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfPileup_PlusCS')

    figure(6)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfPileup_CIE')

    figure(7)
    bar(bins(1:end-1),CntsOfArrived(1:end-1),'histc')
    title('CntsOfPileup_PlusCS_CIE')

    figure(8)
    hold off
    plot(bins,CntsOfArrived,'DisplayName','ED')
    hold on
    plot(bins,CntsOfCIE,'DisplayName','ED + CIE')
    plot(bins,CntsOfPlusCS_CSonly,'DisplayName','CS events + CIE')
    plot(bins,CntsOfPileup,'DisplayName','Pileup')
    plot(bins,CntsOfPileup_PlusCS_PlusCIE,'DisplayName','All effects')
    %plot(bins,CntsOfPileup_CIE,'DisplayName','CntsOfPileup_CIE')
%    plot(bins,CntsOfPileup_PlusCS_CIE,'DisplayName','CntsOfPileup_PlusCS_CIE')
    hold off

    legend

figure(9)
    hold off
    plot(bins,CntsOfArrived/max(CntsOfArrived(2:end)),'-k','LineWidth',2,'DisplayName','Energy Depositions (ED)')
    hold on
    plot(bins,CntsOfCIE/max(CntsOfCIE(2:end)),'-r','LineWidth',2,'DisplayName','ED + CIE')
    plot(bins,CntsOfPlusCS_CSonly/max(CntsOfPlusCS_CSonly(2:end)),'-b','LineWidth',2,'DisplayName','CS events + CIE')
    plot(bins,CntsOfPileup/max(CntsOfPileup(2:end)),'-g','LineWidth',2,'DisplayName','ED + Pileup')
    plot(bins,CntsOfPileup_PlusCS_PlusCIE/max(CntsOfPileup_PlusCS_PlusCIE(2:end)),':m','LineWidth',2,'DisplayName','ED + all effects')
    %plot(bins,CntsOfPileup_CIE/max(CntsOfPileup_CIE),'DisplayName','CntsOfPileup_CIE')
    %plot(bins,CntsOfPileup_PlusCS_CIE/max(CntsOfPileup_PlusCS_CIE),'DisplayName','CntsOfPileup_PlusCS_CIE')
    hold off

    legend
   xlabel('Signal intentisy')
    ylabel('Counts')
xlim([ss, 35])
    
    

figure(10)
    hold off
    plot(bins,CntsOfArrived,'-k','LineWidth',2,'DisplayName','Energy Depositions (ED)')
    hold on
    plot(bins,CntsOfCIE,'--r','LineWidth',2,'DisplayName','ED + CIE')
    plot(bins,CntsOfPlusCS_CIE,'--b','LineWidth',2,'DisplayName','ED + CS + CIE')
    plot(bins,CntsOfPileup,'--g','LineWidth',2,'DisplayName','ED + Pileup')
    plot(bins,CntsOfPileup_PlusCS_PlusCIE,':k','LineWidth',2,'DisplayName','ED + all effects')
    hold off

    legend

    xlabel('Signal intentisy')
    ylabel('Counts')
xlim([18, 35])



    figure(11)
    hold off
    plot(bins,CntsOfArrived/max(CntsOfArrived(18/ss:35/ss)),'-k','LineWidth',2,'DisplayName','Energy Depositions (ED')
    hold on
    plot(bins,CntsOfCIE/max(CntsOfCIE(18/ss:35/ss)),'--r','LineWidth',2,'DisplayName','ED + CIE')
    plot(bins,CntsOfPlusCS_CIE/max(CntsOfPlusCS_CIE(18/ss:35/ss)),'--b','LineWidth',2,'DisplayName','ED + CS + CIE')
    plot(bins,CntsOfPileup/max(CntsOfPileup(18/ss:35/ss)),'--g','LineWidth',2,'DisplayName','ED + Pileup')
    plot(bins,CntsOfPileup_PlusCS_PlusCIE/max(CntsOfPileup_PlusCS_PlusCIE(18/ss:35/ss)),':m','LineWidth',2,'DisplayName','ED + all effects')
    hold off

xlim([18, 35])
    xlabel('Signal intentisy')
    ylabel('Counts')

    legend


    figure(12)
    subplot(3, 3, 2);
    hold off
    plot(bins,CntsOfArrived/max(CntsOfArrived(2:35/ss)),'-k','LineWidth',2,'DisplayName','Energy Depositions (ED')
    hold on
xlim([ss, 35])
    xlabel('Signal intentisy')
    ylabel('Counts')
    
    subplot(3, 3, 4);
    hold off
    plot(bins,CntsOfCIE/max(CntsOfCIE(2:35/ss)),'-r','LineWidth',2,'DisplayName','ED + CIE')
    hold on
xlim([ss, 35])
    xlabel('Signal intentisy')
    ylabel('Counts')

        subplot(3, 3, 5);
    hold off
    plot(bins,CntsOfPileup/max(CntsOfPileup(2:35/ss)),'-g','LineWidth',2,'DisplayName','ED + Pileup')
    hold on
xlim([ss, 35])
    xlabel('Signal intentisy')
    ylabel('Counts')

        subplot(3, 3, 6);
    hold off
    plot(bins,CntsOfPlusCS_CSonly/max(CntsOfPlusCS_CSonly(2:35/ss)),'-b','LineWidth',2,'DisplayName','CS events + CIE')
    hold on
xlim([ss, 35])
    xlabel('Signal intentisy')
    ylabel('Counts')

        subplot(3, 3, 8);
    hold off
    plot(bins,CntsOfPileup_PlusCS_PlusCIE/max(CntsOfPileup_PlusCS_PlusCIE(2:35/ss)),':m','LineWidth',2,'DisplayName','ED + all effects')
    hold on
xlim([ss, 35])
    xlabel('Signal intentisy')
    ylabel('Counts')



 figure(13)
 t = tiledlayout(3,3);

% Add a big shared label for the whole figure
xlabel(t, 'Ouput signal','FontSize',18,'FontWeight','bold');
ylabel(t, 'Counts','FontSize',18,'FontWeight','bold');
%title(t, 'Overall Title for the Figure');

  nexttile(2);
    plot(bins,CntsOfArrived,'-k','LineWidth',2,'DisplayName','Deposited energy spectrum')
    xlim([ss, 35])
    %legend('Location','NorthWest')
    title('Deposited energy spectrum')

  nexttile(4);
    plot(bins,CntsOfCIE,'-r','LineWidth',2,'DisplayName','Induced in-pixel signals')
    xlim([ss, 35])
    %legend
    title('Induced in-pixel signals')

nexttile(5);
    plot(bins,CntsOfPileup,'-g','LineWidth',2,'DisplayName','Pileup effects')
    xlim([ss, 35])
    %legend
    title('Pileup effects')

  nexttile(6);
    plot(bins,CntsOfPlusCS_CSonly,'-b','LineWidth',2,'DisplayName','Induced neighbour-pixel signals')
    xlim([ss, 35])
    %legend
    title('Induced neighbour-pixel signals')

  nexttile(8);
    plot(bins,CntsOfArrived/max(CntsOfArrived(2:end)),'-k','LineWidth',2,'DisplayName','Input spectrum')
    hold on
    plot(bins,CntsOfPileup_PlusCS_PlusCIE/max(CntsOfPileup_PlusCS_PlusCIE(2:end)),':m','LineWidth',2,'DisplayName','Output spectrum')
    xlim([ss, 35])
    %legend('Location','NorthWest')
    title('Normalised input and output spectra')
    


 figure(14)
 t = tiledlayout(2,2);

% Add a big shared label for the whole figure
xlabel(t, 'Ouput signal','FontSize',32,'FontWeight','bold');
ylabel(t, 'Counts','FontSize',32,'FontWeight','bold');
%title(t, 'Overall Title for the Figure');


  nexttile(2);
    plot(bins,CntsOfCIE,'-r','LineWidth',3,'DisplayName','Induced in-pixel signals')
    xlim([ss, 35])
    legend
    title('Induced in-pixel signals','FontSize',24,'FontWeight','bold')
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

nexttile(3);
    plot(bins,CntsOfPileup,'-g','LineWidth',3,'DisplayName','Pileup effects')
    xlim([ss, 35])
    legend
    title('Pileup effects','FontSize',24,'FontWeight','bold')
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

  nexttile(4);
    plot(bins,CntsOfPlusCS_CSonly,'-b','LineWidth',3,'DisplayName','Induced neighbour-pixel signals')
    xlim([ss, 35])
    legend
    title('Induced neighbour-pixel signals','FontSize',24,'FontWeight','bold')
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

  nexttile(1);
    plot(bins,CntsOfArrived/max(CntsOfArrived(2:end)),'-k','LineWidth',3,'DisplayName','Input spectrum')
    hold on
    plot(bins,CntsOfPileup_PlusCS_PlusCIE/max(CntsOfPileup_PlusCS_PlusCIE(2:end)),':m','LineWidth',3,'DisplayName','Output spectrum')
    xlim([ss, 35])
    legend('Location','NorthWest')
    title('Normalised input and output spectra','FontSize',24,'FontWeight','bold')
set(gca, 'FontSize', 14, 'FontWeight', 'bold');


    Raw = sum(CntsOfArrived)
    CS = sum(CntsOfPlusCS_CSonly)
    CIE = sum(CntsOfCIE)
    Pileup = sum(CntsOfPileup)

  %%  
end
end

