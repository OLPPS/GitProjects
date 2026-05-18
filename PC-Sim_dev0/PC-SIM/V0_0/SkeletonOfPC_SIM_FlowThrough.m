%MC Data formatting and set standard file names

%Choose if doing CIE, basic charge cloud projection or lookup table for
%responses (is this different from CIE? perhaps need LUT options? Perhaps use something to generate a LUT?)

%CIE Map formatting and set standard file names

%Generate formatting info and selection info: perhaps including a 
%walkthrough GUIDE Script, like the scintillator type, for PC-SIM to
%read in and format itself with.

%Then--- techinically PC-SIM starts here?
    %INITIALISATION OF SCANNER BUILD
        %Import detector
        %Read in time formatting information
        %Splice MC data to extract relevant timeslice data (looping over all time slices)
     %SIGNAL GENERATION STEPS
            %Perform pixelisation
            %Perform dimensional scaling (if option needed)
            %Generate PseudoEvents (CIE or other method, if option needed)
            %Perform LUT (may involve projection ones too)
    %SIGNAL PROCESSING

    "Moved sinograms folder into main section and added template script for the overall flow with PC-SIM and this folder plan."