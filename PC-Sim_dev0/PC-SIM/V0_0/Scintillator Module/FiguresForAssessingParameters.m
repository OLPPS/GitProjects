%Summary of input parameters and impact on Scintillator-SiPM outputs
clear
clc


for ij = 1:100

    load("ScintCrystalParams.mat")
    load("SiPMParams.mat")

%Scintillator decay
excitationEnergy = 100; %in keV
NumberOfPhotons = round(sum( max(0, Exci2Opti_mean + Exci2Opti_std .* randn(uint16(excitationEnergy), 1))) * SurvivalFraction * PDE );
EmmissionTimes = 1 + exprnd(1/ExponentialDecayConstants(1), NumberOfPhotons, 1);

figure(1)
hold off
histogram(EmmissionTimes)
title("Scintillator decay profile (100 keV photon)")
xlabel('time (ns)')
ylabel('Number of optical photons')

%SPAD recharge
SingleSPADTraceLength = 100;
SPADCharge = zeros(SingleSPADTraceLength,1);
SPADCharge(2:end) = exp(-(1/ExponentialDecayConstants_SPADsRecharge).*[0:SingleSPADTraceLength-2]);
figure(2)
hold off
plot(SPADCharge)
title("SPAD recharge behaviour")
xlabel('time (ns)')
ylabel('SPAD discharge level')

%Offset SPAD recharges (show differences in gain)
%Calculate variations in gain
SPADsGainMap = SPADsGainNoise_Mean + SPADsGainNoise_STD .* randn(SPADsInX,SPADsInY);
if CapVariation == "Y"
lowercap = max(SPADsGainNoise_Mean - SPADsGainNoise_STD * SPADsGainNoise_Cap,0);
uppercap = SPADsGainNoise_Mean + SPADsGainNoise_STD * SPADsGainNoise_Cap;
SPADsGainMap = min(max(SPADsGainMap,lowercap),uppercap);
else
SPADsGainMap = max(SPADsGainMap,0);
end

NumberOfSPADsToSample = 100;
xes = ceil(SPADsInX * rand(NumberOfSPADsToSample,1));
yes = ceil(SPADsInY * rand(NumberOfSPADsToSample,1));

figure(3)
hold off

for iijj = 1:NumberOfSPADsToSample
SingleSPADTraceLength = 100;
SPADCharge = zeros(SingleSPADTraceLength,1);
SPADCharge(2:end) = SPADsGainMap(xes(iijj),yes(iijj)) .* exp(-(1/ExponentialDecayConstants_SPADsRecharge).*[0:SingleSPADTraceLength-2]);
plot(SPADCharge)
title("SPAD recharge behaviour with different gains")
hold on
end
hold off
xlabel('time (ns)')
ylabel('SPAD gain (normalised to mean gain)')



%Multi-Discharge in single SPAD (show capped peaking)
%%%%% THIS IS NOT OPTIMISED FOR LARGE SCALE DEPLOYMENT
SingleSPADTraceLength = 100+1;
NumberOfTriggers = 5;
SPADCharge = zeros(SingleSPADTraceLength,1);
TriggerTimes = ceil(SingleSPADTraceLength * rand(NumberOfTriggers,1))
TriggerTimes = sort(TriggerTimes); 
for iijj = 1:NumberOfTriggers
SPADCharge(1+TriggerTimes(iijj):SingleSPADTraceLength) = exp(-(1/ExponentialDecayConstants_SPADsRecharge).*[0:size([1+TriggerTimes(iijj):SingleSPADTraceLength],2)-1])';
end
figure(4)
hold off
plot(SPADCharge)
title("SPAD recharge behaviour, multiple events")
xlim([0 SingleSPADTraceLength])
ylim([0 1.25])
xlabel('time (ns)')
ylabel('SPAD discharge level')

%Integrated signal across triggers
TotalSignal = zeros(SingleSPADTraceLength,1);
%logicalIncreaseMask = SPADCharge(1:end-1) < SPADCharge(2:end);
%shifts = find(logicalIncreaseMask) + 1;
%signals = SPADCharge(shifts) - SPADCharge(logicalIncreaseMask);
for iijj = 1:NumberOfTriggers
    if TriggerTimes(iijj) < length(SPADCharge)
signals(iijj) = SPADCharge(1+TriggerTimes(iijj))-SPADCharge(TriggerTimes(iijj));
    else
signals(iijj) = 0;
    end
TotalSignal(TriggerTimes(iijj):end) = TotalSignal(TriggerTimes(iijj):end) + signals(iijj);
end
figure(5)
hold off
plot(TotalSignal)
title("Integrated signal from multiple firings")
xlabel('time (ns)')
ylabel('Integral of SPAD discharges')
xlim([0 100])

%Multi-SPADs triggered in single SiPM
X = rand(NumberOfPhotons,1);
Y = rand(NumberOfPhotons,1);
GlobalTrace = SPADsTraceGlobal(X,Y,SPADsInX,SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge);
GlobalTrace_DiffGains = SPADsTraceGlobal_DifferentialGains(X,Y,SPADsInX,SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge,SPADsGainMap);

figure(6)
hold off
plot(GlobalTrace)
title("Multi-SPAD discharge in single SiPM from 100 keV photon")
hold on
plot(GlobalTrace_DiffGains,'--r')
legend('Uniform response','Disperse gains')
hold off
xlabel('time (ns)')
ylabel('Signal from SiPM')


%Combination of all effects
figure(7)
histogram(EmmissionTimes,"BinWidth",1,"Normalization","countdensity")
hold on
plot(GlobalTrace,'-r','LineWidth',1.5)
title("Time profiles of optical photons vs signal out for 100 keV photon")
legend('Optical Photons','Electronic signals')
hold off
xlabel('time (ns)')
ylabel('Optical photons, Signal from SiPM')
ylim([0 250])



GlobalTrace2500 = SPADsTraceGlobal(X,Y,SPADsInX,SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge);

SPADsInX = 10;
SPADsInY = 100;
GlobalTrace1000 = SPADsTraceGlobal(X,Y,SPADsInX,SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge);

SPADsInX = 20;
SPADsInY = 200;
GlobalTrace4000 = SPADsTraceGlobal(X,Y,SPADsInX,SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge);

%SPADsInX = 100;
%SPADsInY = 100;
%GlobalTrace10000 = SPADsTraceGlobal(X,Y,SPADsInX,SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge);

figure(8)
hold off
histogram(EmmissionTimes,"BinWidth",1,"Normalization","countdensity")
title("Time profiles of optical photons vs signal out for  different SPAD densities")
hold on
plot(GlobalTrace1000,'-r','LineWidth',1.5)
plot(GlobalTrace2500,'-k','LineWidth',1.5)
plot(GlobalTrace4000,'-g','LineWidth',1.5)
%plot(GlobalTrace10000,'-c','LineWidth',1.5)
legend('Photons','1000 SPADS','2500 SPADs', '4000 SPADs')
hold off
xlabel('time (ns)')
ylabel('Optical photons, Signal from SiPM')
ylim([0 250])







pause(0.5)
end

%%


%For each event
 %   For each spad
  %      Generate 

%%  NOTE THAT WITH LATER PHOTONS -> LESS SPAD PARTIAL RECHARGING SO MORE ACCURATE COMPARED TO LOSES IN SIGNAL AT HIGHER PHOTON DENSITIES... CAN VARY NUMBER OF SPADS TO DEMONSTRATE THIS (50X5 VS 100X100)

%%  ALSO NOTE ISSUE AT LINE 80 WHICH OCCASIONALLY CAUSES AN OUT OF RANGE VALUE




function GlobalTrace = SPADsTraceGlobal(X,Y,SPADsInX, SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge)
%Build Pixelated Array of SPADs
SPADsArray = accumarray([ceil(X*SPADsInX) ceil(Y*SPADsInY)], sort(EmmissionTimes), [SPADsInX SPADsInY], @(v){v});
GlobalTrace = zeros(max(ceil(EmmissionTimes))+1,1);
for xx = 1:SPADsInX
    for yy = 1:SPADsInY
    TriggerTimes = ceil(SPADsArray{xx,yy});
    signals = zeros(length(TriggerTimes),1);

    SingleSPADTraceLength = max(ceil(EmmissionTimes))+1;
    NumberOfTriggers = length(TriggerTimes);
    SPADCharge = zeros(SingleSPADTraceLength,1);

for iijj = 1:NumberOfTriggers
SPADCharge(1+TriggerTimes(iijj):SingleSPADTraceLength-1) = exp(-(1/ExponentialDecayConstants_SPADsRecharge).*[0:size([1+TriggerTimes(iijj):SingleSPADTraceLength],2)-2])';
end



TotalSignal = zeros(SingleSPADTraceLength,1);
for iijj = 1:NumberOfTriggers
signals(iijj) = SPADCharge(1+TriggerTimes(iijj))-SPADCharge(TriggerTimes(iijj));
TotalSignal(TriggerTimes(iijj):end) = TotalSignal(TriggerTimes(iijj):end) + signals(iijj);
end

GlobalTrace(TriggerTimes) = GlobalTrace(TriggerTimes) + signals;
    end
end

end







function GlobalTrace = SPADsTraceGlobal_DifferentialGains(X,Y,SPADsInX, SPADsInY,NumberOfPhotons,EmmissionTimes,ExponentialDecayConstants_SPADsRecharge,SPADsGainMap)
%Build Pixelated Array of SPADs
SPADsArray = accumarray([ceil(X*SPADsInX) ceil(Y*SPADsInY)], sort(EmmissionTimes), [SPADsInX SPADsInY], @(v){v});
GlobalTrace = zeros(max(ceil(EmmissionTimes))+1,1);
for xx = 1:SPADsInX
    for yy = 1:SPADsInY
    SPADsGain =  SPADsGainMap(xx,yy);
    TriggerTimes = ceil(SPADsArray{xx,yy});
    signals = zeros(length(TriggerTimes),1);

    SingleSPADTraceLength = max(ceil(EmmissionTimes))+1;
    NumberOfTriggers = length(TriggerTimes);
    SPADCharge = zeros(SingleSPADTraceLength,1);

for iijj = 1:NumberOfTriggers
SPADCharge(1+TriggerTimes(iijj):SingleSPADTraceLength-1) = exp(-(1/ExponentialDecayConstants_SPADsRecharge).*[0:size([1+TriggerTimes(iijj):SingleSPADTraceLength],2)-2])';
end



TotalSignal = zeros(SingleSPADTraceLength,1);
for iijj = 1:NumberOfTriggers
signals(iijj) = SPADsGain * (SPADCharge(1+TriggerTimes(iijj))-SPADCharge(TriggerTimes(iijj)));
TotalSignal(TriggerTimes(iijj):end) = TotalSignal(TriggerTimes(iijj):end) + signals(iijj);
end

GlobalTrace(TriggerTimes) = GlobalTrace(TriggerTimes) + signals;
    end
end

end