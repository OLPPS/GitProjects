DEBUGGING = 0;
UseDIPH = 0;
AUTOREMOVER = 1;
saveClassic = 0;
StreamlinedOutput = 1;
TIMING_RESOLUTION = 0;
non_impulsePulseShape = 2;
SAVEPIXBOARD = 1;

SourceDetectorOrientation = 'rev';

CIEIntegrationTime = 100;
FDR = 1;
PCS =1;
DLR =1; 
SR = 1;
IDEAL = 1;
DLR_INTEGRATIONTIME = 10*10^-9;
SR_INTEGRATIONTIME = 10*10^-9;
PCS_INTEGRATIONTIME = 10*10^-9;
FDR_RESETTIME = 90*10^-9;
FDR_RESETDURATION = 10*10^-9;
Add3x3Dy = 1;
Sub3x3Dy = 1;
Add3x3St = 1;
Sub3x3St = 1;
AddHybridDy = 1;
SubHybridDy = 1;
AddHybridSt = 1;
SubHybridSt = 1;
Add2x2Dy = 1;
Sub2x2Dy = 1;
Add2x2St = 1;
Sub2x2St = 1;

DLR_WriteTime = 0;
SR_WriteTime = 0;
    CSCA_SearchingTime = 10*10^(-9);

Sub2x2St_Counters = 0;
Add2x2St_Counters = 0;
Sub2x2Dy_Counters = 0;
Add2x2Dy_Counters = 0;
Sub3x3St_Counters = 0;
Add3x3St_Counters = 0;
Sub3x3Dy_Counters = 0;
Add3x3Dy_Counters = 0;
SubHybridSt_Counters = 0;
AddHybridSt_Counters = 0;
SubHybridDy_Counters = 0;
AddHybridDy_Counters = 0;
NoCSCA_STD_Counters = 0;
NoCSCA_PCS_Counters = 0;
NoCSCA_DLR_Counters = 0;
NoCSCA_SR_Counters = 0;
NoCSCA_FDR_Counters = 0;
NoCSCA_IDEAL_Counters = 0;


    ShapingTime = 1*10^-9;
    TimeWindow = 1*10^-9;
    BoxingTime = 10^(-9);
    
%Temp
mu = 1;
Vb = 1;

THRESHOLDS = [10 30 48 62]; %29 36 44 51 58 64 76];
Thresholds = THRESHOLDS;
NumOfThresh = max(size(Thresholds));




save('PARAMS')