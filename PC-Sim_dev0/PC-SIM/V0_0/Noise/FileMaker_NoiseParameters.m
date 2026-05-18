clear
DeadPixelsFraction = 0.01;
HotPixelFraction = 0.01;
MinimumNoisiness_Healthy = 1;
MinimumNoisiness_Hot = 1;
MaximumNoisiness_Healthy = 4;
MaximumNoisiness_Hot = 4;

NoiseFreq = 10*10^-9;
NoiseFreq2 = NoiseFreq;

outputfilename = "NoiseParameters";

save(outputfilename)