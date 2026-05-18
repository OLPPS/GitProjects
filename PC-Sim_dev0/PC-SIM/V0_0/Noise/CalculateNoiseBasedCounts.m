function ArrayOfCountsToAdd = CalculateNoiseBasedCounts(EventTimesList, Sigma, SigmaTolerance, THRESHOLDS,FREQ2,SimulatedTime)

binStart = EventTimesList(1);
TimesSampled = 1;
dt = FREQ2; %Calculate bin width
        for i = 2:length(EventTimesList)
            if EventTimesList(i) > binStart + dt %DOUBLE CHECK THIS LINE IN WHOLE CONTEXT OF CODE: SHOULD IT BE EventTimesList?
                binStart = EventTimesList(i);  % start a new bin
                TimesSampled = TimesSampled + 1;
            end
        end

%Initialise variable to hold Thresholded data
fractionAboveThreshold = zeros(length(THRESHOLDS),1);
ArrayOfCountsToAdd = zeros(length(THRESHOLDS),1);

%%Determine the number of times noise should be expected to have risen over
%%thershold (within a given tolerance)

        %Calculate number of above threshold random noises expected
        M = 0;
        for T = 1:length(THRESHOLDS)

            if SigmaTolerance >= 0
                fractionAboveThreshold(T) = max(normcdf(SigmaTolerance*Sigma,M,Sigma)- normcdf(THRESHOLDS(T),M,Sigma),0); %%REVISIT THIS
            else
                fractionAboveThreshold(T) = 1- normcdf(THRESHOLDS(T),M,Sigma);

            end
            NumberOfExpectedSamplesAboveThreshold(T) = fractionAboveThreshold(T) * SimulatedTime/FREQ2;
        end

        
        

        for T = 1:length(THRESHOLDS)
            %%Calculate number of events that need to be added
            EventsToBeAdded = NumberOfExpectedSamplesAboveThreshold(T) - (TimesSampled*fractionAboveThreshold(T));

            %%Add events to output array
            ArrayOfCountsToAdd(T) = round(EventsToBeAdded);
        end


end
