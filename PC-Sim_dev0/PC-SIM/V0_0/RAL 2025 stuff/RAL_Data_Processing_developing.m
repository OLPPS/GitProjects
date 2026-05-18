%% To setup datafile access
clear
clc


%{
% Create a MAT-file to store the large array
outfile = 'RAL_data_LF45keV.mat';
DataFile = matfile(outfile, 'Writable', true);
%}

%% To create zeros arrays
%{
% Preallocate the array in the MAT-file
DataFile.data = zeros(1024, 80, 80, 51, 2,'uint32');
%}

%% To overwrite and array with zeros
%{
infile = 'RAL_data_zeros.mat';
DataFile = matfile(infile, 'Writable', true);
for iijj = 1:51
    iijj
    DataFile.data(:,:,:,:,iijj) = zeros(1024,80,80,51,1,'uint32');
end
%}

%%
%{
for iijj = 0:50
    iijj
if iijj < 10
    RawDataFile = strcat('sensor4_45keV_600v_noAttn_20C_19sHisto_50umAS_final_184546_00000',num2str(iijj),('.h5'))
else
    RawDataFile = strcat('sensor4_45keV_600v_noAttn_20C_19sHisto_50umAS_final_184546_0000',num2str(iijj),('.h5'))
end

DataFile.data(:,:,:,:,iijj+1) = h5read(RawDataFile,'/dummy');

end
finished = 1;
%}
%%
clc
clear %EnergyFrame energySpectrum_AllPixInRange tempFrames tempframe

%{
lowery = 46;
uppery = 64;
lowerx = 65;
upperx = 74;
numOfFile = 50;
scanPositionsWithinFile = 50;
noiseFloor = 1;
infile = 'RAL_data_LF45keV.mat';
%}
load("Variables_HF15keV.mat")

DataFile = matfile(infile, 'Writable', false);
for iijj = 1:1:numOfFile+1 %files

    EnergyFrame(:,:,:,:)= DataFile.data(:,lowery:uppery,lowerx:upperx,:,iijj);

energySpectrum_AllPixInRange(:,:) = sum(sum(EnergyFrame,3),2);


tempFrames = sum(EnergyFrame(noiseFloor:end,:,:,:),1);
for ij = 1:1:scanPositionsWithinFile+1 %Position within scan
    %{
        figure(1)
plot([200:1024],energySpectrum_AllPixInRange(200:end,ij),'-o')
hold on
LL = ;
UL = 466;
plot([LL,LL],[0,3000],'-')
plot([UL,UL],[0,3000],'-')
hold off
ylim([0,3100]);
movegui('west')
title(strcat('RowNum = ', num2str(iijj), ' ColNum = ', num2str(ij)))
    %}
figure(7)
tempframe(:,:) = tempFrames(1,:,:,ij); %tempFrames(1,:,:,50) - tempFrames(1,:,:,1);
surf(tempframe)
view(0,90)
title(strcat('RowNum = ', num2str(iijj), ' ColNum = ', num2str(ij)))
%pause(0.1)
movegui('east')
end
end


%446 466
