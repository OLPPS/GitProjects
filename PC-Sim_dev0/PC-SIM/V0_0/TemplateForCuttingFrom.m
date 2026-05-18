ForFemke_IFPb_FDRRedo_ScriptControlled2(0,0,'10e6','FullPhantom_Smaller_RealisticContrast',50)

function ForFemke_IFPb_FDRRedo_ScriptControlled2(lower_bound, upper_bound, IDSTRING, BASEFILENAME, SPLITS)
if nargin == 0
    args = argv;
    lower_bound = str2double(args{1});
    upper_bound = str2double(args{2});
    IDSTRING = args{3};
    BASEFILENAME = args{4};
    SPLITS = str2double(args{5});
    myscript(lower_bound, upper_bound, IDSTRING, BASEFILENAME, SPLITS);
end
%% Skeleton plan
%% Reset workspace
%clear
%clc
%cd /home/mmiteam/CERN24/data_analysis/data_HF250_1_5mm/output
ProcessStage = 0;
GA = datetime;
%% Set up the table of times for this RUN
STARTTIME = 0;
TOTALTIME = 1; %Total time of simulation RUN ON THIS SERVER
PROJECTIONS = 2000;%Total number of projections INCLUDED IN THIS DATASET (this server)

PROJECTIONSPERFILE = PROJECTIONS / SPLITS;
SPLITTIMES(1:SPLITS,1:PROJECTIONSPERFILE+1) = 0; %Array contains start indexes for
%each projection plus last time included in the file (boundaries of each
%split adjacent to each other).
SPLITINDEXES(1:SPLITS,1:PROJECTIONSPERFILE+1) = 0;

STRIPPEDGATEFILENAMES(1:SPLITS) = "";

for RUN_ID = lower_bound:upper_bound%SPLITS - 1 %Set this to the start and end of the batch being processed here (5 scripts in parallel on each of current servers, 2024).


STRIPPEDGATEFILENAMES(RUN_ID+1) = strcat(BASEFILENAME,'_',num2str(RUN_ID),'_',IDSTRING,'.hits.txt_stripped');


for nj = 1:PROJECTIONSPERFILE+1
        splittime = STARTTIME + TOTALTIME/PROJECTIONS * ((RUN_ID * PROJECTIONSPERFILE)+(nj-1));
    SPLITTIMES(RUN_ID+1,nj) = splittime;%find(Eventsdata(:, 4) >= splittime, 1);
end

end







for RUN_ID = lower_bound:upper_bound%SPLITS - 1







SPLITINDEXES(RUN_ID+1,1) = 1;
for ds = 1:PROJECTIONSPERFILE-1
  SPLITINDEXES(RUN_ID+1,ds+1) = find(Eventsdata(:, 4) >= SPLITTIMES(RUN_ID+1,ds+1), 1);
end
SPLITINDEXES(RUN_ID+1,PROJECTIONSPERFILE+1) = size(Eventsdata,1);




%clear Eventsdata
%%


for ProjectionFlag = 1:PROJECTIONSPERFILE
    clearvars -except ProjectionFlag GateDataFile firstprojection secondprojection GEOMETRIES RUN_ID      DEBUGGING UseDIPH CIEIntegrationTime AUTOREMOVER saveClassic TIMING_RESOLUTION FDR non_impulsePulseShape PCS DLR SR SourceDetectorOrientation DLR_INTEGRATIONTIME SR_INTEGRATIONTIME PCS_INTEGRATIONTIME FDR_RESETTIME FDR_RESETDURATION Add3x3Dy Sub3x3Dy Add3x3St Sub3x3St AddHybridDy SubHybridDy AddHybridSt SubHybridSt Add2x2Dy Sub2x2Dy Add2x2St Sub2x2St xlength ylength zlength AvX AvY AvZ Xcom Ycom Zcom CIEMapsAllTimes RUNSTART RUNSTOP RUNDURATION IDEAL Eventsdata SPLITS TOTALTIME PROJECTIONS PROJECTIONSPERFILE SPLITINDEXES BASEFILENAME IDSTRING STRIPPEDGATEFILENAMES SPLITTIMES CIEMapHeader   xlength ylength zlength AvX AvY AvZ Xcom Ycom Zcom CIEMapsAllTimes K GA GB StreamlinedOutput SAVEPIXBOARD 


if NOISEON == 0
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
A = datetime;
%Store xpixnum, ypixnum, time, signal, pileupsignal
%function SPCS_pixelateBoard_v0.0(file,Extension,VERBOSE,DEPTH,pitch,FLUXX,ShapingTime)
%end
%% Apply ACTS
%Define Thresholds
algstart = datetime;
Thresholds = THRESHOLDS;
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
                            if PixelatedBoard(X,Y,1,pointer+1) - PixelatedBoard(X,Y,1,pointer) < 10^-9
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













elseif NOISEON == 1

    
    
    
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
A = datetime;
%Store xpixnum, ypixnum, time, signal, pileupsignal
%function SPCS_pixelateBoard_v0.0(file,Extension,VERBOSE,DEPTH,pitch,FLUXX,ShapingTime)
%end
%% Apply ACTS
%Define Thresholds
algstart = datetime;
Thresholds = THRESHOLDS;
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
                            if PixelatedBoard(X,Y,1,pointer+1) - PixelatedBoard(X,Y,1,pointer) < 10^-9
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


algstop = datetime;
algtime = algstop - algstart

%%
        strcat("Saving projection number ", num2str(ProjectionFlag), " of ", num2str(PROJECTIONSPERFILE),".")
        output_label = strcat(LABEL,'_projectionData');
        save(output_label, 'GEOMETRIES', 'Sub2x2St_Counters', 'Add2x2St_Counters', 'Sub2x2Dy_Counters', 'Add2x2Dy_Counters', 'Sub3x3St_Counters', 'Add3x3St_Counters', 'Sub3x3Dy_Counters', 'Add3x3Dy_Counters', 'SubHybridSt_Counters', 'AddHybridSt_Counters', 'SubHybridDy_Counters', 'AddHybridDy_Counters', 'NoCSCA_STD_Counters', 'NoCSCA_PCS_Counters', 'NoCSCA_DLR_Counters', 'NoCSCA_SR_Counters', 'NoCSCA_FDR_Counters', 'NoCSCA_IDEAL_Counters', 'RUN_ID', 'DEBUGGING', 'UseDIPH', 'CIEIntegrationTime', 'AUTOREMOVER', 'saveClassic', 'StreamlinedOutput', 'BASEFILENAME', 'IDSTRING', 'TIMING_RESOLUTION', 'FDR', 'non_impulsePulseShape', 'PCS', 'DLR', 'SR', 'SourceDetectorOrientation', 'DLR_INTEGRATIONTIME', 'SR_INTEGRATIONTIME', 'PCS_INTEGRATIONTIME', 'FDR_RESETTIME', 'FDR_RESETDURATION', 'Add3x3Dy', 'Sub3x3Dy', 'Add3x3St', 'Sub3x3St', 'AddHybridDy', 'SubHybridDy', 'AddHybridSt', 'SubHybridSt', 'Add2x2Dy', 'Sub2x2Dy', 'Add2x2St', 'Sub2x2St', 'IDEAL', 'DECVALS', 'FWTM','-v7.3')

clearvars -except IDEAL ProcessStage SPLITS TOTALTIME PROJECTIONS PROJECTIONSPERFILE SPLITINDEXES BASEFILENAME IDSTRING STRIPPEDGATEFILENAMES       LABEL Sub2x2St_Counters Add2x2St_Counters Sub2x2Dy_Counters Add2x2Dy_Counters Sub3x3St_Counters Add3x3St_Counters Sub3x3Dy_Counters Add3x3Dy_Counters SubHybridSt_Counters AddHybridSt_Counters SubHybridDy_Counters AddHybridDy_Counters NoCSCA_STD_Counters NoCSCA_PCS_Counters NoCSCA_DLR_Counters NoCSCA_SR_Counters NoCSCA_FDR_Counters NoCSCA_IDEAL_Counters GEOMETRIES Eventsdata ProjectionFlag GEOMETRIES RUN_ID      DEBUGGING UseDIPH CIEIntegrationTime AUTOREMOVER saveClassic TIMING_RESOLUTION FDR non_impulsePulseShape PCS DLR SR SourceDetectorOrientation DLR_INTEGRATIONTIME SR_INTEGRATIONTIME PCS_INTEGRATIONTIME FDR_RESETTIME FDR_RESETDURATION Add3x3Dy Sub3x3Dy Add3x3St Sub3x3St AddHybridDy SubHybridDy AddHybridSt SubHybridSt Add2x2Dy Sub2x2Dy Add2x2St Sub2x2St SPLITTIMES CIEMapHeader      xlength ylength zlength AvX AvY AvZ Xcom Ycom Zcom CIEMapsAllTimes K GA GB StreamlinedOutput SAVEPIXBOARD FDR_RESETTIME
output_label = strcat(LABEL,'_processed')

if StreamlinedOutput == 1
else
    save(output_label,'-v7.3')
end

GB = datetime;

C = GB - GA

end
end
end


end
%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%

%%
%Checked orientation of GATE vs COMSOL (positive spatial diretion is towards cathode, vs positive spatial direction is towards anode). Correction applied to GATE data so that now positive spatial direction is towards anode. Thus, when the negative shift is produced by the addition of a positive offset PRIOR to the inversion of coordinates, the result is that all points of interaction are moved towards the cathode direction. Therefore, photons which penetrate further into the pixel on GATE (anode end, higher energy) may remain in the pixel after the shift, whilst the lower energy photons (those nearer the cathode end initially) are shifted outside of the pixel. When these are excluded, this thus appears as a beam hardening effect.
