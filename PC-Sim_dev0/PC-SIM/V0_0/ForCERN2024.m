%% Skeleton plan
%% Reset workspace
clear
clc
ProcessStage = 0;
%% Setup error flags and other handler parameters
for ErrorFlags = 1:1
    PixilationInXErrorExit = 0;
    PixilationInYErrorExit = 0;
    PixilationInZErrorExit = 0;
    OverwriteX = 0;
    OverwriteY = 0;
    OverwriteZ = 0;
    CatFail = 0;
    SkippedRuns = 0;
    ERROR = 0;
    IndicesChecker(1:9) = 0;
    CIEEqualedNaN = 0;
    CIEZeros = 0;
    XExtrapolationFlag = 0;
    YExtrapolationFlag = 0;
    ZExtrapolationFlag = 0;
end

  

%% Input parameters (should be read in eventually)
for InputParameters = 1:1
    ShapingTime = 1*10^-9
    TimeWindow = 1*10^-9
    VERBOSE = 1;
    if VERBOSE == 1
        fprintf('Reading in geometry and physics process parameters...');
        startA = datetime;
        for ReadInGeometries = 1:1
            PixelWidth = GEOMETRIES(1); %mm
            PixelPitch = GEOMETRIES(2); %mm
            PixelDepth = GEOMETRIES(3); %mm
            xOffset = GEOMETRIES(4); %mm
            yOffset = GEOMETRIES(5); %mm
            zOffset = GEOMETRIES(6); %mm %NOT YET TESTED AS NOT YET NEEDED
            MagXEdgeSpace = GEOMETRIES(7); %mm %NOT YET TESTED AS NOT YET NEEDED
            MagYEdgeSpace = GEOMETRIES(8); %mm %NOT YET TESTED AS NOT YET NEEDED
            PixelsInXDirection = GEOMETRIES(9); %pixels - 1
            PixelsInYDirection = GEOMETRIES(10); %pixels - 1
        end

        for ReadInPhysicsProcesses = 1:1
            ChargeSharing = PhysicsProcesses(1);
            FiniteCCSize = PhysicsProcesses(2);
            Detrapping = PhysicsProcesses(3);
            e_eRepulsion = PhysicsProcesses(4);
            ElectronicNoise = PhysicsProcesses(5);
            non_impulsePulseShape = PhysicsProcesses(6);
            DefectivePixels = PhysicsProcesses(7);
        end
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

end

for FreeingWorkspaceMemory = 1:1
%Free up some memory
    if VERBOSE == 1
        fprintf('Freeing workspace...')
        startA = datetime;
    end

    clear ErrorFlags

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

end

%% Take in COMSOL map
%Temp
mu = 1;
Vb = 1;
DecayConstant = 4385000;
THRESHOLDS = [0 10 50 100];
BaseName = 'CdTe_100_30_075_1000V_3mm-1'
for i = 1:9
    fprintf('CIEMapAddress %i...', i);
%{  
    switch i
        case 1
            LocalisedAddress = '_centre_actual';
        case 2
            LocalisedAddress = '_centre_adjacent';
        case 3
            LocalisedAddress = '_centre_diagonal';
        case 4
            LocalisedAddress = '_edge_actual';
        case 5
            LocalisedAddress = '_edge_adjacent';
        case 6
            LocalisedAddress = '_edge_diagonal';
        case 7
            LocalisedAddress = '_corner_actual';
        case 8
            LocalisedAddress = '_corner_adjacent';
        case 9
            LocalisedAddress = '_corner_diagonal';
    end
            
        MapAddress{i} = strcat(BaseName,LocalisedAddress);
        CIEMapAddressesFORNOW{i} = ProcessZero(MapAddress{i});
fprintf('Done\n')
%}
        CIEMapAddresses{i} = '1mmCdTe100um30ns_MERGED.txt';
end
%Temp
for FreeingWorkspaceMemory = 1:1
%Free up some memory
    if VERBOSE == 1
        fprintf('Freeing workspace...')
        startA = datetime;
    end

    clear InputParameters

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

end


%% Take in any charge sharing maps

%% Take in any other maps


%% Take in raw GATE data then remove out of pixel events, 0 energy events and convert to intrapixel coordinates
for ReadInGATEData = 1:1
    if VERBOSE == 1
        fprintf('Reading in Gate data...')
        startA = datetime;
    end
    GateDataFile = dlmread(RawDataAddress);

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end
end

for FormattingGATEData = 1:1
    if VERBOSE ==1
        fprintf('Identifying zero energy events and outside strikes from Gate data...')
        startA = datetime;
    end
    for CountingAndSizingGateData = 1:1
        NumberOfEvents = size(GateDataFile,1);
        Indices(1:NumberOfEvents) = 0;
        energyzeroes = 0;
        OutsideStrikes = 0;
        for ColumnLabels = 1:1
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
        end


        for W = 1:NumberOfEvents
            if GateDataFile(W,energy) == 0
                energyzeroes = energyzeroes +1;
            else
                xpostemp = GateDataFile(W,xpos) + xOffset - MagXEdgeSpace;
                ypostemp = GateDataFile(W,ypos) + yOffset - MagYEdgeSpace;
                if xpostemp < 0 || ypostemp < 0 || xpostemp > (PixelPitch * (PixelsInXDirection+1)) || ypostemp > (PixelPitch * (PixelsInYDirection+1))
                    OutsideStrikes = OutsideStrikes + 1;
                else
                    Indices(W) = W;
                end
            end
        end
        NumberOfEventsRaw = NumberOfEvents;
        NumberOfEvents = NumberOfEvents - energyzeroes - OutsideStrikes;
        EventsOnlyOutputData = coder.nullcopy(zeros(NumberOfEvents,Columns));

    end

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

    if VERBOSE == 1
        fprintf('Constructing Event data array of non-zero energy events within range...')
        startA = datetime;
    end

    for FilterEnergyZeroAndOutstrikeEventsAndPixelate = 1:1
        jj = 0;
        JustEventIndices = nonzeros(Indices);
        for M = 1:size(JustEventIndices,1)
            W = JustEventIndices(M);
            jj = jj + 1;
            xpostemp = GateDataFile(W,xpos) + xOffset - MagXEdgeSpace;
            ypostemp = GateDataFile(W,ypos) + yOffset - MagYEdgeSpace;
            %zpos = text_fileA(W,zpos) + PixelDepth/2 + zOffset; Hard
            %coded in zpos assignation

            %Calculate pixel number
            EventsOnlyOutputData(jj,xpixnum) = fix(xpostemp/PixelPitch);
            EventsOnlyOutputData(jj,ypixnum) = fix(ypostemp/PixelPitch);

            %Normalise to distance within pixel
            EventsOnlyOutputData(jj,xpos) = (xpostemp - (round(EventsOnlyOutputData(jj,xpixnum))*PixelPitch))/PixelPitch;
            EventsOnlyOutputData(jj,ypos) = (ypostemp - (round(EventsOnlyOutputData(jj,ypixnum))*PixelPitch))/PixelPitch;
            EventsOnlyOutputData(jj,zpos) = 1.0 - ((GateDataFile(W,zpos) + PixelDepth/2 + zOffset)/ PixelDepth); %This is to account for z being reversed between Gate and Comsol
            EventsOnlyOutputData(jj,time) = GateDataFile(W,time);
            EventsOnlyOutputData(jj,energy) = (GateDataFile(W,energy));

        end


    end

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

    for ReshapingEventsArrayForUse = 1:1
           if VERBOSE == 1
        fprintf('Correct for edge rounding errors and sort pixels by type...')
        startA = datetime;
    end 

    %Reshape to account for xpos = 0 but xpixnum = 210 (rounding
        %error)
        for W = 1:max(size(EventsOnlyOutputData))
            if EventsOnlyOutputData(W,xpixnum) == PixelsInXDirection +1
                EventsOnlyOutputData(W,xpixnum) = EventsOnlyOutputData(W,xpixnum) - 1;
                EventsOnlyOutputData(W,xpos) = 1 - EventsOnlyOutputData(W,xpos);
            end
            if EventsOnlyOutputData(W,ypixnum) == PixelsInYDirection +1
                EventsOnlyOutputData(W,ypixnum) = EventsOnlyOutputData(W,ypixnum) - 1;
                EventsOnlyOutputData(W,ypos) = 1 - EventsOnlyOutputData(W,ypos);
            end
        end

        %%Classify pixels as central, edge or corner type
        for W = 1:max(size(EventsOnlyOutputData))
            %if x = 0 or x = 255 or y = 0 or y = 255 in a 256 x 256 pixel array
            if (round(EventsOnlyOutputData(W,xpixnum)) == 0) || (round(EventsOnlyOutputData(W,xpixnum)) == PixelsInXDirection)
                EventsOnlyOutputData(W,pixtype) = EventsOnlyOutputData(W,pixtype) + 1;
            end
            if (round(EventsOnlyOutputData(W,ypixnum)) == 0) || (round(EventsOnlyOutputData(W,ypixnum)) == PixelsInYDirection)
                EventsOnlyOutputData(W,pixtype) = EventsOnlyOutputData(W,pixtype) + 1;
            end
        end
    
    
    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end
    
    end
end

for FreeingWorkspaceMemory = 1:1
%Free up some memory
    if VERBOSE == 1
        fprintf('Freeing workspace...')
        startA = datetime;
    end

    clear GateDataFile Indices JustEventIndices M W xpostemp ypostemp FreeingWorkspaceMemory FormattingGATEData ReadInGATEData CountingAndSizingGateData FilterEnergyZeroAndOutstrikeEventsAndPixelate ReshapingEventsArrayForUse

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

end

if DEBUGGING == 1
fprintf('Saving Output Data For Debug Logs...');
startA = datetime;
save(strcat(LABEL,'_DBG1'), '-v7.3')
stopB = datetime;
fprintf('Completed in (hr:min:sec): %s\n', stopB-startA);
end

%% DIPH module could possibly go here


%% Sort events based on which CIE maps would be needed
for SortingEventsByCIEMapToUse = 1:1
    if VERBOSE ==1
        fprintf('Sorting events based on CIE map needed...')
        startA = datetime;
    end

    %Sort Gate Data Into 9 Separate Lists (CAc,CAd,CDi, EAc, EAd, EDi, CoAc,
    %CoAd and CoDi)
    Lists(1:9,1:NumberOfEvents) = 0; %THIS MAY BE AN ISSUE WHICH CAUSES A MEMORY CRASH
    Indices(1:9) = 0;
    Map2Loc(1:9) = 0;
    %Count how many events exist in each list type
    for k = 1:NumberOfEventsToCalculate
        CIEMapTypeToUse = 3*EventsOnlyOutputData(k,pixtype)+(EventsOnlyOutputData(k, DriftState)+1); %DriftState 0 = event, 1 = adjacent, 2 = diagonal
        Map2Loc(CIEMapTypeToUse) = Map2Loc(CIEMapTypeToUse) + 1;
        Lists(CIEMapTypeToUse,Map2Loc(CIEMapTypeToUse)) = k;
    end
    Indices(1:9) = Map2Loc(1:9); %set list length for each CIE map type
    IndicesChecker(1:9) = Indices(1:9); %backup list lengths

    %build each list with ID numbers for each event
    ListOfCentresActual(1:Indices(1)) = Lists(1,1:Indices(1));
    ListOfCentresAdjacent(1:Indices(2)) = Lists(2,1:Indices(2));
    ListOfCentresDiagonal(1:Indices(3)) = Lists(3,1:Indices(3));
    ListOfEdgesActual(1:Indices(4)) = Lists(4,1:Indices(4));
    ListOfEdgesAdjacent(1:Indices(5)) = Lists(5,1:Indices(5));
    ListOfEdgesDiagonal(1:Indices(6)) = Lists(6,1:Indices(6));
    ListOfCornersActual(1:Indices(7)) = Lists(7,1:Indices(7));
    ListOfCornersAdjacent(1:Indices(8)) = Lists(8,1:Indices(8));
    ListOfCornersDiagonal(1:Indices(9)) = Lists(9,1:Indices(9));

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA);
    end

end

%% Calculate CIE values
for CalculateCIEValues = 1:1

    if VERBOSE == 1
        startA = datetime;
        fprinft('Calculating CIE values for...')
    end    
    
%Only consider maps that are likely to be needed.
%%NOTE: THIS BIT COULD BE REMOVED IF THERE IS NO DIPH USED
if UseDIPH == 1
    Lim = 1:9;
else
    Lim = [1, 4, 7];
end

%Calculate CIE values for each list considered in turn
for CaseNumber = Lim
    startC = datetime;
    
    for ErrorCheckAndReset = 1:1
        %Check if errors in pixilation are known and if so skip, otherwise proceed
            if  PixilationInXErrorExit == 1 ||  PixilationInYErrorExit == 1 ||  PixilationInZErrorExit == 1
                SkippedRuns =  SkippedRuns +1;
                %                        save(LABELLA,'-struct', 'StructArray','-v7.3');
                fprintf('Skipped %i runs\n', SkippedRuns)
            else
                clear vals
            end
    end
    
    for IdentifyCorrectList = 1:1
        switch CaseNumber
            case 1
                vals = ListOfCentresActual;
                PixelType = 'ListOfCentresActual';
            case 2
                vals = ListOfCentresAdjacent;
                PixelType = 'ListOfCentresAdjacent';
            case 3
                vals = ListOfCentresDiagonal;
                PixelType = 'ListOfCentresDiagonal';
            case 4
                vals = ListOfEdgesActual;
                PixelType = 'ListOfEdgesActual';
            case 5
                vals = ListOfEdgesAdjacent;
                PixelType = 'ListOfEdgesAdjacent';
            case 6
                vals = ListOfEdgesDiagonal;
                PixelType = 'ListOfEdgesDiagonal';
            case 7
                vals = ListOfCornersActual;
                PixelType = 'ListOfCornersActual';
            case 8
                vals = ListOfCornersAdjacent;
                PixelType = 'ListOfCornersAdjacent';
            case 9
                vals = ListOfCornersDiagonal;
                PixelType = 'ListOfCornersDiagonal';
        end
    end
    
    disp("\nPixels of type... " + PixelType)
    
    %Import and build correct Map
    [xlength, ylength, zlength, AvX, AvY, AvZ, Xcom, Ycom, Zcom, CIEMapsAllTimes] = BuildAllCIEMapsNeeded(CIEMapHeader, CIEIntegrationTime);
 
    K(1:xlength,1:ylength,1:zlength) = 0;
    for ii = 1:xlength
        for jj = 1:ylength
            for kk = 1:zlength

               K(ii,jj,kk) = max(CIEMapsAllTimes(ii,jj,kk,:)); %SHOULD THIS BE MINIMUM FOR ADJACENT PIXELS WITH NEGATIVE SIGNAL INDUCTION? CONSIDER, BUT PROBABLY NOT AS CHARGE SHARING LIKELY TO DOMINATE? PERHAPS THIS SHOULD BE A SEPARATE MODULE? SHOULD THIS BE ON NEGATIVE FOR ADJACENT AND DIAGONAL PIXELS WHEN WE HAVE CHARGE CLOUD SIZE CLACULATIONS VIA MAGI, TO ALLOW FOR APPLICATIONS LIKE HEXITEC?

            end
        end
    end


    %Go through relevant part of Events, based on relevant Lists, and
    %calculate CIEs
    for CalculateCIEOfEventsInCurrentLists = 1:1
        %Move through points and interpolate IF RELEVANT TO THIS MAP
        x1 = 0;
        x2 = 0;
        y1 = 0;
        y2 = 0;
        z1 = 0;
        z2 = 0;
        for i = vals
            tx =  EventsOnlyOutputData(i, zpos)*10^4; %Swapping x for z to convert from Gate to Comsol coordiante systems
            ty =  EventsOnlyOutputData(i, ypos)*10^4; %The *10^4 here is hard coded assuming the AvX, AvY and AvZ values used an integer variable with the same scaling for the mapping key
            tz =  EventsOnlyOutputData(i, xpos)*10^4; %Swapping z for x to convert from Gate to Comsol coordiante systems


            %The 10^4 in this section is hard coded assuming the int64 format key is
            %being used
            for LocateCornersOfCIEMapVolumeContainingTheEventBeingAssessed = 1:1 %function [x1,x2,y1,y2,z1,z2] = Locator(StepX,StepY,StepZ, targetx, targety, targetz, lengthx,lengthy,lengthz)
                LoopContinue = 0;
                %Similar for each spatial dimension variable
                x1 = fix(tx/AvX)+1; %Calculate lower bound
                x2 = x1+1; %Increment to upper bound
                if x1 > xlength-1 %make any necessary adjustments if events over upper edge
                    if x1 > xlength + 1 %UNTESTED CAPTURE FLAG...
                        ZExtrapolationFlag = ZExtrapolationFlag+1;
                        EventsOnlyOutputData(i, CIE) = -41572;
                        LoopContinue = 1;
                    end %... TO HERE
                    x1 = xlength - 1;
                    x2 = xlength;
                elseif x1 < 1 %make any necessary adjustments if events below lower edge
                    if x1 < 0 %UNTESTED CAPTURE FLAG...
                        ZExtrapolationFlag = ZExtrapolationFlag+1;
                        EventsOnlyOutputData(i, CIE) = -41572;
                        LoopContinue = 1;
                    end %... TO HERE
                    x1 = 1;
                    x2 = 2;
                end

                y1 = fix(ty/AvY)+1; %Calculate lower bound
                y2 = y1+1; %Increment to upper bound
                if y1 > ylength-1
                    if y1 > ylength + 1
                        YExtrapolationFlag = YExtrapolationFlag+1;
                        EventsOnlyOutputData(i, CIE) = -41572;
                        LoopContinue = 1;
                    end
                    y1 = ylength - 1;
                    y2 = ylength;
                elseif y1 < 1
                    if y1 < 0
                        YExtrapolationFlag = YExtrapolationFlag+1;
                        EventsOnlyOutputData(i, CIE) = -41572;
                        LoopContinue = 1;
                    end
                    y1 = 1;
                    y2 = 2;
                end

                z1 = fix(tz/AvZ)+1; %Calculate lower bound
                z2 = z1+1; %Increment to upper bound
                if z1 > zlength -1
                    if z1 > zlength +1
                        XExtrapolationFlag = XExtrapolationFlag+1;
                        EventsOnlyOutputData(i, CIE) = -41572;
                        LoopContinue = 1;
                    end
                    z1 = zlength - 1;
                    z2 = zlength;
                elseif z1 < 1
                    if z1 < 0
                        XExtrapolationFlag = XExtrapolationFlag+1;
                        EventsOnlyOutputData(i, CIE) = -41572;
                        LoopContinue = 1;
                    end
                    z1 = 1;
                    z2 = 2;
                end

            end

            for TrilinearInterpolation = 1:1
                if LoopContinue == 1;
                    continue
                end

                %Locate CIE values for the 8 points surrounding the event being assessed
                P111 = K(x1,y1,z1);
                P112 = K(x1,y1,z2);
                P121 = K(x1,y2,z1);
                P122 = K(x1,y2,z2);
                P211 = K(x2,y1,z1);
                P212 = K(x2,y1,z2);
                P221 = K(x2,y2,z1);
                P222 = K(x2,y2,z2);

                %Extract the spatial coordinates of the 8 points surrounging the event
                %being assessed
                x1val = double((Xcom(x1)));
                x2val = double((Xcom(x2)));
                y1val = double((Ycom(y1)));
                y2val = double((Ycom(y2)));
                z1val = double((Zcom(z1)));
                z2val = double((Zcom(z2)));

                %Interpolate from a cube to a square
                ratios(1:6) = 0;
                ratios(1) = (x2val-tx)/(x2val-x1val);
                ratios(2) = (tx-x1val)/(x2val-x1val);
                TP1 = (P111*ratios(1))+(P211*ratios(2));
                TP2 = (P112*ratios(1))+(P212*ratios(2));
                TP3 = (P121*ratios(1))+(P221*ratios(2));
                TP4 = (P122*ratios(1))+(P222*ratios(2));

                %Interpolate from a square to a line
                ratios(3) = (y2val-ty)/(y2val-y1val);
                ratios(4) = (ty-y1val)/(y2val-y1val);
                TP5 = (TP1*ratios(3))+(TP3*ratios(4));
                TP6 = (TP2*ratios(3))+(TP4*ratios(4));

                %Interpolate from a line to a point, and assign the event that CIE value
                EventsOnlyOutputData(i, CIE) = (TP5*(z2val-tz)/(z2val-z1val))+(TP6*(tz-z1val)/(z2val-z1val));

                %Catch for CIE not being equal to a number. This could indicate a
                %serious error and may need to be investigated.
                if isnan( EventsOnlyOutputData(i, CIE)) == 1
                    fprintf('Error in calculating CIEs: CIE = NaN\nSetting CIE value to zero\n')
                    fprintf('%i, %i, %i, %i, %i, %i, %i \n', x1val, x2val,y1val,y2val,z1val,z2val, tz);
                    pause(10)
                    CIEEqualedNaN = CIEEqualedNaN + 1;
                    EventsOnlyOutputData(i, CIE) = 0;
                end
                if EventsOnlyOutputData(i,CIE) == 0
                    CIEZeros = CIEZeros + 1;
                end

            end

        end



    end
    stopD = datetime;
    fprintf('...Completed in (hr:min:sec): %s\n', stopD-startA);

    if VERBOSE == 1
        stopB = daetime;
        fprintf('\nCompleted in (hr:min:sec): %s\n', stopB-startA);
    end
end

end

for FreeingWorkspaceMemory = 1:1
%Free up some memory
    if VERBOSE == 1
        fprintf('Freeing workspace...')
        startA = datetime;
    end

    clear CaseNumber startC stopD ErrorCheckAndReset IdentifyCorrectList CalculateCIEValues 
    clear IdentifyCorrectList vals PixelType CalculateCIEOfEventsInCurrentLists
    clear tx ty tz x1 x2 y1 y2 z1 z2 LocateCornersOfCIEMapVolumeContainingTheEventBeingAssessed
    clear P111 P112 P121 P122 P211 P212 P221 P222 x1val x2val y1val y2val z1val z2val TP1 TP2 TP3
    clear TP4 TP5 TP6 ratios TrilinearInterpolation
    clear ListOfCentresActual ListOfCentresAdjacent ListOfCentresDiagonal ListOfEdgesActual ListOfEdgesAdjacent
    clear ListOfEdgesDiagonal ListOfCornersActual ListOfCornersAdjacent ListOfCornersDiagonal
    clear Lists Indices Event Map2Loc
    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

end

%% Remove 0 CIE values
for RemoveCIEZeros = 1:1
    NumberOfEvents = size(EventsOnlyOutputData,1);
    Indices(1:NumberOfEvents) = 0;
    CIEzeroes = 0;
    for W = 1:NumberOfEvents
        if  EventsOnlyOutputData(W,CIE) == 0
            CIEzeroes = CIEzeroes +1;
        else
            Indices(W) = W;
        end
    end

    if AUTOREMOVER == 1
        RemoveZeroCIEs = AUTOREMOVER;
    else
        fprintf('%i events detected with CIE = 0.\n This represents approximately %i percent of the data.\n Would you like to remove them?\n', CIEZeros, 100*CIEZeros/(NumberOfEvents))
        RemoveZeroCIEs = input('Yes = 1\n No = 0\n')
        %%NOTE: SHOULD PUT IN HERE SOMETHING ABOUT EXITING IF ANSWER IS NOT
        %%A 1 OR A 0
    end
    if RemoveZeroCIEs == 1

        %Remove events that do not deposit any energy
        fprintf('\nRemoving any events with CIE = 0\n')
        fprintf('Copying Data for processing... \n')

        AllEventsOnlyOutputData = EventsOnlyOutputData;
        clear EventsOnlyOutputData
        NumberOfEventsRaw2 = NumberOfEvents;
        NumberOfEvents = NumberOfEvents - CIEzeroes;
        fprintf('Done.\nCreating Data array for holding results...\n')
        EventsOnlyOutputData = coder.nullcopy(zeros(NumberOfEvents,Columns));
        JustEventIndices = nonzeros(Indices);
        clear Indices
        fprintf('Done.\nConstructing new array with relevant events...\n')

        Sz = size(JustEventIndices,1);

        EventsOnlyOutputData(1:Sz,:) = AllEventsOnlyOutputData(JustEventIndices(:),:);
    else
        Sz = NumberOfEvents;
    end

end

for FreeingWorkspaceMemory = 1:1
%Free up some memory
    if VERBOSE == 1
        fprintf('Freeing workspace...')
        startA = datetime;
    end

    clear AllEventsOnlyOutputData CIEzeros AllEventsOnlyOutputData

    if VERBOSE == 1
        stopB = datetime;
        fprintf('Completed in (hr:min:sec): %s\n', stopB-startA)
    end

end

%% Optionally save as classic CoGI only EventsOnlyOutput data output format
    if saveClassic == 1
        save(strcat(LABEL,'_ClassicCoGIAndFormat'), '-v7.3')
    end

%% Adjust times based on transit times
%THIS IS AN INTERESTING QUESTION. May be irrelvant on the time scales
%considered though?
for CorrectTimesForFiniteDriftTime = 1:1
if DriftTimeCorrections == 1 %% will need to resort events once pixelated, or in EventsOnlyArray or "by hand" as each event is put in to a sorted array. Sort is marginally quicker on sorted than unsorted data
    for ii = 1:Sz
        if CrystalOrientation == 1 
            DriftTime = EventsOnlyOutputData(ii,(1-zpos))/(mu*Vb); % check the 1-zpos. Should it be xpos? or 1-xpos? or zpos etc.? Make a document to help show orientation%(d*l)/(mu*Vb) but l = 1 as normalised spatial coordinates are used.
        else
            DriftTime = EventsOnlyOutputData(ii,(zpos))/(mu*Vb); % check the 1-zpos. Should it be xpos? or 1-xpos? or zpos etc.? Make a document to help show orientation%(d*l)/(mu*Vb) but l = 1 as normalised spatial coordinates are used.
        end
        EventsOnlyOutputData(ii,time) = EventsOnlyOutputData(ii,time) + DriftTime;
    end
end
end

if DEBUGGING == 1
fprintf('Saving Output Data For Debug Logs...');
startA = datetime;
save(strcat(LABEL,'_DBG2'), '-v7.3')
stopB = datetime;
fprintf('Completed in (hr:min:sec): %s\n', stopB-startA);
end

%% pixelatedBoard(CIE adjusted events in a pixelated grid)
%For building, debugging and testing, allow to be called either as a
%file address or data array (EventsOnly list form) directly.

XMax = PixelsInXDirection;
YMax = PixelsInYDirection;

EventsOnlyOutputData(:,xpixnum) = EventsOnlyOutputData(:,xpixnum) +1;
EventsOnlyOutputData(:,ypixnum) = EventsOnlyOutputData(:,ypixnum) +1;

% Calculate Bin sizes needed
Counter(1:XMax,1:YMax) = 0; 
for N = 1 : Sz
    Counter(EventsOnlyOutputData(N,xpixnum),EventsOnlyOutputData(N,ypixnum)) = Counter(EventsOnlyOutputData(N,xpixnum),EventsOnlyOutputData(N,ypixnum)) +1;
end
MaxBinDepth = max(max(Counter))+2;

%%Build PixelatedBoard and other tools
for BuildingPixelatedBoardAndTools = 1:1

PixelatedBoard(1:XMax,1:YMax,1:3,1:MaxBinDepth) = 0;
PixelatedBoard(1:XMax,1:YMax,1,1:MaxBinDepth) = 2.0*max(EventsOnlyOutputData(:,time));
Headers(1:XMax,1:YMax,1:2) = 2; %This ensures first value has time = 0; May need to account for this if doing a pixelated N ~= MasterList N test later? Equally means there will be NumOfPixels more events than expected if using headers(2) to calcualte events without subtracting 1 from headers(2) first.
MasterList(1:MaxSize,1:4) = -1000;

%PixelatedBoard
Tyme = 1;
ProtoSig = 2;
HashNum = 3;
PrevSignal = 4;

%Headers
Pointer = 1;
ListSize = 2;

%MasterList
NNum = 1;
Xval = 2;
Yval = 3;
Tag = 4;
end

%%Populate PixelatedBoard
for PopulatePixelatedBoard = 1:1
    for N = 1:MaxSize
        X = EventsOnlyOutputData(N,xpixnum);
        Y = EventsOnlyOutputData(N,ypixnum);
        PixelatedBoard(X,Y,1,Headers(X,Y,2)) = EventsOnlyOutputData(N,time);
        PixelatedBoard(X,Y,2,Headers(X,Y,2)) = EventsOnlyOutputData(N,CIE)*EventsOnlyOutputData(N,energy)*1000;
        PixelatedBoard(X,Y,3,Headers(X,Y,2)) = N;
        %   PixelatedBoard(X,Y,4,Headers(X,Y,2)) = 1;
        Headers(X,Y,2) = Headers(X,Y,2)+1;
        MasterList(N,1) = N; % Could be removed?
        MasterList(N,2) = X;
        MasterList(N,3) = Y;
        MasterList(N,4) = 1;
    end

        clear EventsOnlyOutputData
    
    MasterList_Backup = MasterList;

end

if DEBUGGING == 1
fprintf('Saving Output Data For Debug Logs...');
startA = datetime;
save(strcat(LABEL,'_DBG3'), 'PixelatedBoard', 'MasterList_Backup', '-v7.3')
stopB = datetime;
fprintf('Completed in (hr:min:sec): %s\n', stopB-startA);
end

%% Decay correct signals to be a running total
here you should do something like the noCSCA case already done, but with no search window (so just each pixel, whole list, using header 2 perhaps? Could still use only a single header) but adding decayed old event to current event.

for X = 1: XMax
    for Y = 1:YMax
        for pointer = 2:Headers(X,Y,2)-1
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

%Store xpixnum, ypixnum, time, signal, pileupsignal
%function SPCS_pixelateBoard_v0.0(file,Extension,VERBOSE,DEPTH,pitch,FLUXX,ShapingTime)
%end
%% Apply CSCAs
%Define Thresholds
Thresholds = THRESHOLDS;
NumOfThresh = max(size(Thresholds));
% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
NoCSCA_Counters(1:XMax, 1:YMax,1:NumOfThresh) = 0;
NoCSCA_EventsProcessed = 0;
for X = 1:XMax
    for Y = 1:YMax
        BelowThresh = 0;
        for pointer = 2:Headers(X,Y,2)-1
            for Thresh = 1:NumOfThresh
                if PixelatedBoard(X,Y,4,pointer) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                    BelowThresh = Thresh;
                    break
                end
            end
            if BelowThresh ~= 0
                AboveThresh = 0;
                for Thresh = BelowThresh:NumberOfThresh
                    if PixelatedBoard(X,Y,2,pointer) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                        AboveThresh = Thresh;
                    else
                        break
                    end
                end
                if AboveThresh ~= 0
                        NoCSCA_Counters(X,Y,BelowThresh:AboveThresh) = NoCSCA_Counters(X,Y,BelowThresh:AboveThresh) + 1;
                end
            end
            NoCSCA_EventsProcessed = NoCSCA_EventsProcessed + 1;
            if rem(NoCSCA_EventsProcessed,Sz/100) == 0
                NoCSCA_PercentageCompleted = NoCSCA_EventsProcessed/(Sz/100)
            end
        end
    end
end

%% 3x3 Dynamic additive
% FOR CSCAs: IN HERE SOMEWHERE YOU NEED TO FIND AN EVENT AND THEN CHECK THE CHARGE ON ALL RELEVANT ADJACENT PIXELS AS OLD HEADER TIME, TO NOW, DECAY CORRECTED, OR SIMILAR.
Add3x3Dy_Counters(1:XMax, 1:YMax,1:NumOfThresh) = 0;
Add3x3Dy_EventsProcessed = 0;
MasterList = MasterList_Backup;
Headers(:,:,1) = 2;

for nn = 1:Sz
    if MasterList(nn,4) == 1
        CSCA_Trigger = 0;
        X = MasterList(2);
        Y = MasterList(3);
        NewLocalSignal = PixelatedBoard(X,Y,2,Headers(X,Y,1));
        OldLocalSignal = PixelatedBoard(X,Y,4,Headers(X,Y,1));
        for Thresh = 1:NumOfThresh
            if NewLocalSignal > Thresholds(Thresh) && OldLocalSignal < Thresholds(Thresh)
                CSCA_Trigger = 1;
                break
            end
        end

        if CSCA_Trigger == 1
        % DO CSCA based THRESHOLD COUNTING
            Signal = NewLocalSignal;
            for dx = -1:1
                for dy = -1:1
                if X+dx >= 1 && X+dx <= XMax && Y+dy >= 1 && Y+dy <= YMax
                    if dx == 0 && dy == 0 
                    else
                        deltaT = PixelatedBoard(X,Y,1,Headers(X,Y,1)) - PixelatedBoard(X+dx,Y+dy,1,Headers(X,Y,1)-1);
                        Signal = Signal + PixelatedBoard(X+dx,Y+dy,2,Headers(X,Y,1)-1)*exp(-DecayConstant*deltaT);
                    end
                end
                end
            end
Need to advance to cover all events within next shaping time?
e.g. 
%Do thresholding on CSCA bit
HOW DO I ADJUST FOR EVENTS BEING RESUMMED?

        end
    end
end



for X = 1:XMax
    for Y = 1:YMax
        BelowThresh = 0;
        for pointer = 2:Headers(X,Y,2)-1
            for Thresh = 1:NumOfThresh
                if PixelatedBoard(X,Y,4,pointer) < Thresholds(Thresh)   %if OLD is below then assign BELOWT = THRESH
                    BelowThresh = Thresh;
                    break
                end
            end
            if BelowThresh ~= 0
                AboveThresh = 0;
                for Thresh = BelowThresh:NumberOfThresh
                    if PixelatedBoard(X,Y,2,pointer) > Thresholds(Thresh) %NEW is above then assign AboveT = Thresh and next
                        AboveThresh = Thresh;
                    else
                        break
                    end
                end
                if AboveThresh ~= 0
                    for Incrementing = BelowThresh:AboveThresh
                        Add3x3Dy_Counters(X,Y,Incrementing) = Add3x3Dy_Counters(X,Y,Incrementing) + 1;
                    end
                end
            end
            Add3x3Dy_EventsProcessed = Add3x3Dy_EventsProcessed + 1;
            if rem(Add3x3Dy_EventsProcessed,Sz/100) == 0
            Add3x3Dy_PercentageCompleted = Add3x3Dy_EventsProcessed/(Sz/100)
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%

function [xlength, ylength, zlength, AvX, AvY, AvZ, Xcom, Ycom, Zcom, CIEMapsAllTimes] = BuildAllCIEMapsNeeded(CIEMapHeader, CIEIntegrationTime)
for formatto4x4AndNormaliseCIEMap =1:1
    PixilationInXErrorExit = 0;
    PixilationInYErrorExit = 0;
    PixilationInZErrorExit = 0;

    CurrentCIEMapAddress = strcat(CIEMapHeader, '1ns_PreProcessed.txt'); % '_', num2str(1), 'ns.txt');

    list = dlmread(CurrentCIEMapAddress);
    list(isnan(list))=0;
    
    rows = max(size(list))/4;
    j = 1;
    Y(1:rows,1:4) = 0;
    for i = 1:rows
        Y(i,1) = list(j);
        Y(i,2) = list(j+1);
        Y(i,3) = list(j+2);
        Y(i,4) = list(j+3);
        j = j + 4;
    end

    

    [minx,minidx] = min(Y(1:rows,1));
    [miny,minidy] = min(Y(1:rows,2));
    [minz,minidz] = min(Y(1:rows,3));
    [maxx,maxidx] = max(Y(1:rows,1));
    [maxy,maxidy] = max(Y(1:rows,2));
    [maxz,maxidz] = max(Y(1:rows,3));
    difx = (maxx-minx);
    dify = (maxy-miny);
    difz = (maxz-minz);

    W(1:max(size(Y)),1:3) = Y(:,1:3);

    for i = 1:rows
        W(i,1) = (Y(i,1) - minx) / difx;
        W(i,2) = (Y(i,2) - miny) / dify;
        W(i,3) = (Y(i,3) - minz) / difz;
    end


    for determineDimensionsOfCIEMaps = 1:1 %reproduces dimensionplotter function
        %            [a, lengthx] = DimensionPlotter(W(:,1));
        %            [b, lengthy] = DimensionPlotter(W(:,2));
        %            [c, lengthz] = DimensionPlotter(W(:,3));



        exit2 = 0;
        lengthx = 1;
        VariableSpace = 0;
        a = 0;
        dimension = 0;
        clear dimension
        dimension = W(:,1);
        dimension(max(size(dimension))+1) = -1;
        while exit2 == 0
            if dimension(lengthx+1) == dimension(lengthx)
                lengthx = lengthx+1;
            elseif dimension(lengthx+1) > dimension(lengthx)
                VariableSpace = VariableSpace+1;
                a(VariableSpace) = dimension(lengthx);
                lengthx = lengthx+1;
            elseif dimension(lengthx+1) < dimension(lengthx)
                VariableSpace = VariableSpace+1;
                a(VariableSpace) = dimension(lengthx);
                exit2 = -1;
            end
        end



        exit2 = 0;
        lengthy = 1;
        VariableSpace = 0;
        b = 0;
        clear dimension
        dimension = W(:,2);
        dimension(max(size(dimension))+1) = -1;
        while exit2 == 0
            if dimension(lengthy+1) == dimension(lengthy)
                lengthy = lengthy+1;
            elseif dimension(lengthy+1) > dimension(lengthy)
                VariableSpace = VariableSpace+1;
                b(VariableSpace) = dimension(lengthy);
                lengthy = lengthy+1;
            elseif dimension(lengthy+1) < dimension(lengthy)
                VariableSpace = VariableSpace+1;
                b(VariableSpace) = dimension(lengthy);
                exit2 = -1;
            end
        end



        exit2 = 0;
        lengthz = 1;
        VariableSpace = 0;
        c = 0;
        clear dimension
        dimension = W(:,3);
        dimension(max(size(dimension))+1) = -1;
        while exit2 == 0
            if dimension(lengthz+1) == dimension(lengthz)
                lengthz = lengthz+1;
            elseif dimension(lengthz+1) > dimension(lengthz)
                VariableSpace = VariableSpace+1;
                c(VariableSpace) = dimension(lengthz);
                lengthz = lengthz+1;
            elseif dimension(lengthz+1) < dimension(lengthz)
                VariableSpace = VariableSpace+1;
                c(VariableSpace) = dimension(lengthz);
                exit2 = -1;
            end
        end

    end

    for pixelateCIEMapsForKey = 1:1
        % compare lengthx and xlength to see if they are the
        % same size. May be able to save a few lines below.

        CIEMatrix = Y(:,4);
        Xcom = int16(a*10^4);
        Ycom = int16(b*10^4);
        Zcom = int16(c*10^4);
        xlength = max(size(Xcom));
        ylength = max(size(Ycom));
        zlength = max(size(Zcom));

        % May then clear a, b and c

    end



    %Generate regular pixilation approximately matching that given in teh
    %imported CIE Map. This is needed so that the step size is averaged over
    %the whole range rather than taken directly from one specific step, to
    %reduce the risk of rounding errors. This was added when the pseudo hash
    %feature was employed, to speed up the CIEMap searching.
    for GenerateRegularPixelation = 1:1
        %These values are not int16 as they are sub unity fractions by definition
        AvX = (1.0/(xlength-1))*10^4; %Calculate step size, assuming points are evenly distributed along length and capped at each end
        AvY = (1.0/(ylength-1))*10^4;
        AvZ = (1.0/(zlength-1))*10^4;
        for i = 1:xlength
            Xapprox(i) = int16(AvX*(i-1)); %Populate a list of points to assess
        end
        for i = 1:ylength
            Yapprox(i) = int16(AvY*(i-1));
        end
        for i = 1:zlength
            Zapprox(i) = int16(AvZ*(i-1));
        end
    end
    fprintf('1\n')


    %Check approximation still gives nearest points
    for CheckRegularisation = 1:1
        CatFail = 0;
        for n = 1:1 %Check approximation still gives nearest points
            IsNotRegular = 0;
            for i = 1:xlength
                if Xcom(i) == Xapprox(i)
                elseif abs(Xcom(i) - Xapprox(i)) < AvX %checks if each artificially pixelation is the nearest to the actual point in the CIE map it si meant to represent
                elseif abs(Xcom(i) - Xapprox(i)) < (2*AvX) %allows some tolerance in pixelation mismatch without exxiting program. 2 means artificial pixel is no more than 1 pixel away from its intended location
                    IsNotRegular = IsNotRegular + 1;
                else
                    IsNotRegular = IsNotRegular + 1;
                    CatFail = 1;
                end
            end

            %Interput if pixel approximation is found to deviate beyond tolerance
            %(X)
            while  CatFail == 1
                Overwrite = input('!!Warning! X-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current critical sub level of 2 voxels. \nPlease consider using a coarser pixilisation in future.\nDo you wish to proceed anyway?\n Yes = 1\n No = 0\n')
                if Overwrite == 1
                    CatFail = 0;
                    OverwriteX = 1;
                    PixilationInXErrorExit = 0;
                elseif Overwrite == 0
                    CatFail = 0;
                    PixilationInXErrorExit = 1;
                    OverwriteX = 0;
                end
            end

            %Output warning that approximations may be inadequate, if the user has
            %decided to go ahead anyway (X)
            if IsNotRegular > 0
                IsNotRegular
                fprintf('!!X-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current sub 1 voxel tolerance. \nPlease consider using a coarser pixilisation.\n')
            end


            IsNotRegular = 0;
            for i = 1:ylength
                if Ycom(i) == Yapprox(i)
                elseif abs(Ycom(i) - Yapprox(i)) < AvY
                elseif abs(Ycom(i) - Yapprox(i)) < (2*AvY)
                    IsNotRegular = IsNotRegular + 1;
                else
                    IsNotRegular = IsNotRegular + 1;
                    CatFail = 1;
                end
            end

            %Interput if pixel approximation is found to deviate beyond tolerance (Y)
            while  CatFail == 1
                Overwrite = input('!!Warning! Y-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current critical sub level of 2 voxels. \nPlease consider using a coarser pixilisation in future.\nDo you wish to proceed anyway?\n Yes = 1\n No = 0\n')
                if Overwrite == 1
                    CatFail = 0;
                    OverwriteY = 1;
                    PixilationInYErrorExit = 0;
                elseif Overwrite == 0
                    CatFail = 0;
                    OverwriteY = 0;
                    PixilationInYErrorExit = 1;
                end
            end

            %Output warning that approximations may be inadequate, if the user has
            %decided to go ahead anyway (Y)
            if IsNotRegular > 0
                IsNotRegular
                fprintf('!!Y-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current sub 1 voxel tolerance. \nPlease consider using a coarser pixilisation.\n')
            end

            IsNotRegular = 0;
            for i = 1:zlength
                if Zcom(i) == Zapprox(i)
                elseif abs(Zcom(i) - Zapprox(i)) < AvZ
                elseif abs(Zcom(i) - Zapprox(i)) < (2*AvZ)
                    IsNotRegular = IsNotRegular + 1;
                else
                    IsNotRegular = IsNotRegular + 1;
                    CatFail = 1;
                end
            end

            %Interput if pixel approximation is found to deviate beyond tolerance (Z)
            while  CatFail == 1
                Overwrite = input('!!Warning! Z-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current critical sub level of 2 voxels. \nPlease consider using a coarser pixilisation in future.\nDo you wish to proceed anyway?\n Yes = 1\n No = 0\n')
                if Overwrite == 1
                    CatFail = 0;
                    OverwriteZ = 1;
                    PixilationInZErrorExit = 0;
                elseif Overwrite == 0
                    CatFail = 0;
                    OverwriteZ = 0;
                    PixilationInZErrorExit = 1;
                end
            end
            if  PixilationInZErrorExit == 1 ||  PixilationInXErrorExit == 1 ||  PixilationInYErrorExit == 1
                Skipped = 1;
            end

            %Output warning that approximations may be inadequate, if the user has
            %decided to go ahead anyway (z)
            if IsNotRegular > 0
                IsNotRegular
                fprintf('!!Z-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current sub 1 voxel tolerance. \nPlease consider using a coarser pixilisation.\n')
            end
        end
    end


    %Build Map
    for BuildCIEMap = 1:1
        p = 1; %Steps through the CIE map in the order of incrementing x, then y, then z: this was the COMSOL format observed
        K(1:xlength,1:ylength,1:zlength) = 0;
        CIEMapsAllTimes(1:xlength,1:ylength,1:zlength,1:CIEIntegrationTime) = 0;
        for k = 1:zlength
            for j = 1:ylength
                for i = 1:xlength
                    K(i,j,k) = CIEMatrix(p);
                    p = p+1;
                end
            end
        end
    end

    CIEMapsAllTimes(1:xlength,1:ylength,1:zlength,1) = K(1:xlength,1:ylength,1:zlength);









    
    for TimePos = 2:CIEIntegrationTime


        CurrentCIEMapAddress = strcat(CIEMapHeader,num2str(CIEIntegrationTime),'ns_PreProcessed.txt'); %'_', num2str(TimePos), 'ns.txt'); 
        
        list = dlmread(CurrentCIEMapAddress);
        rows = max(size(list))/4;
        j = 1;
        Y(1:rows,1:4) = 0;
        for i = 1:rows
            Y(i,1) = list(j);
            Y(i,2) = list(j+1);
            Y(i,3) = list(j+2);
            Y(i,4) = list(j+3);
            j = j + 4;
        end

        [minx,minidx] = min(Y(1:rows,1));
        [miny,minidy] = min(Y(1:rows,2));
        [minz,minidz] = min(Y(1:rows,3));
        [maxx,maxidx] = max(Y(1:rows,1));
        [maxy,maxidy] = max(Y(1:rows,2));
        [maxz,maxidz] = max(Y(1:rows,3));
        difx = (maxx-minx);
        dify = (maxy-miny);
        difz = (maxz-minz);

        W(1:max(size(Y)),1:3) = Y(:,1:3);

        for i = 1:rows
            W(i,1) = (Y(i,1) - minx) / difx;
            W(i,2) = (Y(i,2) - miny) / dify;
            W(i,3) = (Y(i,3) - minz) / difz;
        end


        for determineDimensionsOfCIEMaps = 1:1 %reproduces dimensionplotter function
            %            [a, lengthx] = DimensionPlotter(W(:,1));
            %            [b, lengthy] = DimensionPlotter(W(:,2));
            %            [c, lengthz] = DimensionPlotter(W(:,3));



            exit2 = 0;
            lengthx = 1;
            VariableSpace = 0;
            a = 0;
            dimension = 0;
            clear dimension
            dimension = W(:,1);
            dimension(max(size(dimension))+1) = -1;
            while exit2 == 0
                if dimension(lengthx+1) == dimension(lengthx)
                    lengthx = lengthx+1;
                elseif dimension(lengthx+1) > dimension(lengthx)
                    VariableSpace = VariableSpace+1;
                    a(VariableSpace) = dimension(lengthx);
                    lengthx = lengthx+1;
                elseif dimension(lengthx+1) < dimension(lengthx)
                    VariableSpace = VariableSpace+1;
                    a(VariableSpace) = dimension(lengthx);
                    exit2 = -1;
                end
            end



            exit2 = 0;
            lengthy = 1;
            VariableSpace = 0;
            b = 0;
            clear dimension
            dimension = W(:,2);
            dimension(max(size(dimension))+1) = -1;
            while exit2 == 0
                if dimension(lengthy+1) == dimension(lengthy)
                    lengthy = lengthy+1;
                elseif dimension(lengthy+1) > dimension(lengthy)
                    VariableSpace = VariableSpace+1;
                    b(VariableSpace) = dimension(lengthy);
                    lengthy = lengthy+1;
                elseif dimension(lengthy+1) < dimension(lengthy)
                    VariableSpace = VariableSpace+1;
                    b(VariableSpace) = dimension(lengthy);
                    exit2 = -1;
                end
            end



            exit2 = 0;
            lengthz = 1;
            VariableSpace = 0;
            c = 0;
            clear dimension
            dimension = W(:,3);
            dimension(max(size(dimension))+1) = -1;
            while exit2 == 0
                if dimension(lengthz+1) == dimension(lengthz)
                    lengthz = lengthz+1;
                elseif dimension(lengthz+1) > dimension(lengthz)
                    VariableSpace = VariableSpace+1;
                    c(VariableSpace) = dimension(lengthz);
                    lengthz = lengthz+1;
                elseif dimension(lengthz+1) < dimension(lengthz)
                    VariableSpace = VariableSpace+1;
                    c(VariableSpace) = dimension(lengthz);
                    exit2 = -1;
                end
            end

        end

        for pixelateCIEMapsForKey = 1:1
            % compare lengthx and xlength to see if they are the
            % same size. May be able to save a few lines below.

            CIEMatrix = Y(:,4);
            Xcom = int16(a*10^4);
            Ycom = int16(b*10^4);
            Zcom = int16(c*10^4);
            xlength = max(size(Xcom));
            ylength = max(size(Ycom));
            zlength = max(size(Zcom));

            % May then clear a, b and c

        end



        %Generate regular pixilation approximately matching that given in teh
        %imported CIE Map. This is needed so that the step size is averaged over
        %the whole range rather than taken directly from one specific step, to
        %reduce the risk of rounding errors. This was added when the pseudo hash
        %feature was employed, to speed up the CIEMap searching.
        for GenerateRegularPixelation = 1:1
            %These values are not int16 as they are sub unity fractions by definition
            AvX = (1.0/(xlength-1))*10^4; %Calculate step size, assuming points are evenly distributed along length and capped at each end
            AvY = (1.0/(ylength-1))*10^4;
            AvZ = (1.0/(zlength-1))*10^4;
            for i = 1:xlength
                Xapprox(i) = int16(AvX*(i-1)); %Populate a list of points to assess
            end
            for i = 1:ylength
                Yapprox(i) = int16(AvY*(i-1));
            end
            for i = 1:zlength
                Zapprox(i) = int16(AvZ*(i-1));
            end
        end
        fprintf('%i \n', TimePos)


        %Check approximation still gives nearest points
        for CheckRegularisation = 1:1
            CatFail = 0;
            for n = 1:1 %Check approximation still gives nearest points
                IsNotRegular = 0;
                for i = 1:xlength
                    if Xcom(i) == Xapprox(i)
                    elseif abs(Xcom(i) - Xapprox(i)) < AvX %checks if each artificially pixelation is the nearest to the actual point in the CIE map it si meant to represent
                    elseif abs(Xcom(i) - Xapprox(i)) < (2*AvX) %allows some tolerance in pixelation mismatch without exxiting program. 2 means artificial pixel is no more than 1 pixel away from its intended location
                        IsNotRegular = IsNotRegular + 1;
                    else
                        IsNotRegular = IsNotRegular + 1;
                        CatFail = 1;
                    end
                end

                %Interput if pixel approximation is found to deviate beyond tolerance
                %(X)
                while  CatFail == 1
                    Overwrite = input('!!Warning! X-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current critical sub level of 2 voxels. \nPlease consider using a coarser pixilisation in future.\nDo you wish to proceed anyway?\n Yes = 1\n No = 0\n')
                    if Overwrite == 1
                        CatFail = 0;
                        OverwriteX = 1;
                        PixilationInXErrorExit = 0;
                    elseif Overwrite == 0
                        CatFail = 0;
                        PixilationInXErrorExit = 1;
                        OverwriteX = 0;
                    end
                end

                %Output warning that approximations may be inadequate, if the user has
                %decided to go ahead anyway (X)
                if IsNotRegular > 0
                    IsNotRegular
                    fprintf('!!X-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current sub 1 voxel tolerance. \nPlease consider using a coarser pixilisation.\n')
                end


                IsNotRegular = 0;
                for i = 1:ylength
                    if Ycom(i) == Yapprox(i)
                    elseif abs(Ycom(i) - Yapprox(i)) < AvY
                    elseif abs(Ycom(i) - Yapprox(i)) < (2*AvY)
                        IsNotRegular = IsNotRegular + 1;
                    else
                        IsNotRegular = IsNotRegular + 1;
                        CatFail = 1;
                    end
                end

                %Interput if pixel approximation is found to deviate beyond tolerance (Y)
                while  CatFail == 1
                    Overwrite = input('!!Warning! Y-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current critical sub level of 2 voxels. \nPlease consider using a coarser pixilisation in future.\nDo you wish to proceed anyway?\n Yes = 1\n No = 0\n')
                    if Overwrite == 1
                        CatFail = 0;
                        OverwriteY = 1;
                        PixilationInYErrorExit = 0;
                    elseif Overwrite == 0
                        CatFail = 0;
                        OverwriteY = 0;
                        PixilationInYErrorExit = 1;
                    end
                end

                %Output warning that approximations may be inadequate, if the user has
                %decided to go ahead anyway (Y)
                if IsNotRegular > 0
                    IsNotRegular
                    fprintf('!!Y-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current sub 1 voxel tolerance. \nPlease consider using a coarser pixilisation.\n')
                end

                IsNotRegular = 0;
                for i = 1:zlength
                    if Zcom(i) == Zapprox(i)
                    elseif abs(Zcom(i) - Zapprox(i)) < AvZ
                    elseif abs(Zcom(i) - Zapprox(i)) < (2*AvZ)
                        IsNotRegular = IsNotRegular + 1;
                    else
                        IsNotRegular = IsNotRegular + 1;
                        CatFail = 1;
                    end
                end

                %Interput if pixel approximation is found to deviate beyond tolerance (Z)
                while  CatFail == 1
                    Overwrite = input('!!Warning! Z-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current critical sub level of 2 voxels. \nPlease consider using a coarser pixilisation in future.\nDo you wish to proceed anyway?\n Yes = 1\n No = 0\n')
                    if Overwrite == 1
                        CatFail = 0;
                        OverwriteZ = 1;
                        PixilationInZErrorExit = 0;
                    elseif Overwrite == 0
                        CatFail = 0;
                        OverwriteZ = 0;
                        PixilationInZErrorExit = 1;
                    end
                end
                if  PixilationInZErrorExit == 1 ||  PixilationInXErrorExit == 1 ||  PixilationInYErrorExit == 1
                    Skipped = 1;
                end

                %Output warning that approximations may be inadequate, if the user has
                %decided to go ahead anyway (z)
                if IsNotRegular > 0
                    IsNotRegular
                    fprintf('!!Z-Pixilation error!! \nSome points may be incorrectly interpolated by greater \nthan the current sub 1 voxel tolerance. \nPlease consider using a coarser pixilisation.\n')
                end
            end
        end


        %Build Map
        for BuildCIEMap = 1:1
            p = 1; %Steps through the CIE map in the order of incrementing x, then y, then z: this was the COMSOL format observed
            K(1:xlength,1:ylength,1:zlength) = 0;
            for k = 1:zlength
                for j = 1:ylength
                    for i = 1:xlength
                        K(i,j,k) = CIEMatrix(p);
                        p = p+1;
                    end
                end
            end
        end

        CIEMapsAllTimes(1:xlength,1:ylength,1:zlength,TimePos) = K(1:xlength,1:ylength,1:zlength);

    end


end
end

%%
%Checked orientation of GATE vs COMSOL (positive spatial diretion is towards cathode, vs positive spatial direction is towards anode). Correction applied to GATE data so that now positive spatial direction is towards anode. Thus, when the negative shift is produced by the addition of a positive offset PRIOR to the inversion of coordinates, the result is that all points of interaction are moved towards the cathode direction. Therefore, photons which penetrate further into the pixel on GATE (anode end, higher energy) may remain in the pixel after the shift, whilst the lower energy photons (those nearer the cathode end initially) are shifted outside of the pixel. When these are excluded, this thus appears as a beam hardening effect.
