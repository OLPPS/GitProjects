%{
clear
clc

A = BuildDetector("NoiseParameters",23,7,12345)
%}

function [outputFilename] = BuildDetector(NoiseParametersFile,PixelsInXDirection,PixelsInYDirection,RandomSeed)
load(NoiseParametersFile)


%% Calculate per pixel noisiness map 
% Generate a map if all pixels were healthy
NoiseLevelMap_Healthy = CalculateGlobalNoiseSigmaDistribution(MaximumNoisiness_Healthy,MinimumNoisiness_Healthy,PixelsInXDirection,PixelsInYDirection,RandomSeed);

% Generate a map if all pixels were hot
NoiseLevelMap_Hot = CalculateGlobalNoiseSigmaDistribution(MaximumNoisiness_Hot,MinimumNoisiness_Hot,PixelsInXDirection,PixelsInYDirection,RandomSeed*17);

% Generate pixel array with distiriubtion of healths
rng(RandomSeed*231,'twister')
PixelHealth = rand(PixelsInXDirection, PixelsInYDirection);

% Use pixel health to determine location of dead pixels (-1), hot pixels (1) and normal pixels (0)
PixelHealth(PixelHealth < DeadPixelsFraction) = -1;
PixelHealth(PixelHealth > (1-HotPixelFraction)) = 1;
PixelHealth(PixelHealth >= DeadPixelsFraction & PixelHealth <= HotPixelFraction ) = 0;

% Adjust map for pixels based on pixel Health value
NoiseLevelMap_Adjusted = NoiseLevelMap_Healthy;
NoiseLevelMap_Adjusted(PixelHealth == -1) = 0;
NoiseLevelMap_Adjusted(PixelHealth == 1) = NoiseLevelMap_Hot(PixelHealth == 1);

outputFilename = strcat("DetectorHlthAndNoisinessMap_RndSd_",num2str(RandomSeed),"_X",num2str(PixelsInXDirection),"_Y",num2str(PixelsInYDirection));
save(outputFilename);
end



function [GlobalNoiseSigmaDistribution] = CalculateGlobalNoiseSigmaDistribution(MaximumNoisiness,MinimumNoisiness,PixelsInXDirection,PixelsInYDirection,RANDOMSEED)
% Calculate global distribution of pixel noisinesses
MeanGlobal = (MaximumNoisiness + MinimumNoisiness)/2;     % Mean should be 0 as electronic noise is modelled as 0 mean, due to calibration
SigmaGlobal = (MaximumNoisiness - MinimumNoisiness)/6;     % Standard deviation should be set such that the 3 sigma value is the acceptable noise level being used (mean plus 3 sigma should be below the noise floor in a properly calibrated system).

rng(RANDOMSEED, 'twister');
GlobalNoiseSigmaDistribution = MeanGlobal + SigmaGlobal .* randn(PixelsInXDirection, PixelsInYDirection);

% Cap noisiness levels at MinimumNoisieness and MaximumNoisieness
GlobalNoiseSigmaDistribution(GlobalNoiseSigmaDistribution < MinimumNoisiness) = MinimumNoisiness;
GlobalNoiseSigmaDistribution(GlobalNoiseSigmaDistribution > MaximumNoisiness) = MaximumNoisiness;
end
