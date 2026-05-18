%% Will assume images are nested in a hierarchy according to Pitch, Flux, Alg, Approach

%Define reference ROIs obtained from an image of reference pitch

clear
clc

iPitch = 0;

%Set initial directory (root node)
BaseDirectory = 'E:\Femke Work Laptop\Reconstructions\Full_recons\FBP\Thesholded';
NPS = zeros(3,7,18,1,101,101);
NPS(:,:,:,:,:,:) = -743286;

for Pitch = 100:50:400
    Pitch
    iPitch = iPitch + 1;
    iFlux = 0;
    for Flux = 6:8
        Flux
        iFlux = iFlux + 1;
        for Alg = 1:18
            for recon = 1:1
                switch recon
                    case 1
                        Recon = 'FBP';
                    case 2
                        Recon = 'PHP';
                    case 3
                        Recon = 'MLM'
                    case 4
                        Recon = 'MSEM';
                end

                %Calculate and navigate to relevant directory - THESE LINES SHOULD BE
                %MODIFIED BASED ON HIERARCHY STRUCTURE.
                for CalcAndNav_Director = 1:1
                    A = num2str(Pitch);
                    B = num2str(Flux);
                    C = num2str(Alg);
                    D = Recon;

                    TargetDirectory = fullfile(BaseDirectory,strcat('Pitch_',A));

                    cd(TargetDirectory)
                end



                %Calculate image name- THESE LINES SHOULD BE MODIFIED BASED ON FILE NAMING
                %CONVENTIONS
                for CalcAndNav_ImageFile = 1:1

                    ImageName = strcat(Recon,'_10e',num2str(Flux),'_approach_',num2str(Alg),'_t1.png');
                end



                %Load and handle image
                for loadimage = 1:1
                    ImageData = LoadAndProcessImage(fullfile(TargetDirectory,ImageName));
                    figure(1)
                    imshow(ImageData,[]);
                end








                %Calculate and display ROIs
                for CalculatingROIs = 1:1

IodineConcentrations_Actual= [35, 30, 25, 20, 15, 10];

CalciumConcentrations_Actual = [300, 260, 220, 180, 140, 100];


                    %Define ROIs for a given reference image.
                    insertsCalcium_mm = [-20.5,6; -26.5,6; -20.5,0; -26.5,0; -20.5, -6; -26.5, -6]/2; % in mm
                    insertsIodine_mm = [20.5,6; 26.5,6; 20.5,0; 26.5,0; 20.5,-6; 26.5,-6]/2; % in mm
                    insertsMixed_mm = [-6,-20.5; -6,-26.5; 0,-20.5; 0,-26.5; 6,-20.5; 6,-26.5]/2; % in mm
                    insertsTungsten_mm = [0, 12.5];

                    pixel_size_1 = Pitch/1000; % mm/pixel
                    width_phantom_mm = 32; % mm
                    %width_phanotm_pixels = 255; % pixels for pitch 0.25
                    width_phantom_pixels = 255/(pixel_size_1/0.25);
                    calculated_pixel_size = width_phantom_mm/width_phantom_pixels;
                    magnification = pixel_size_1/calculated_pixel_size;

                    centre_p = size(ImageData, 1)/2 + 0.5;
                    centresCalcium = insertsCalcium_mm/calculated_pixel_size + centre_p;
                    centresIodine = insertsIodine_mm/calculated_pixel_size + centre_p;
                    centresMixed = insertsMixed_mm/calculated_pixel_size + centre_p;
                    centresTungsten = insertsTungsten_mm/calculated_pixel_size +centre_p;

                    sizes = 2*magnification/(2*pixel_size_1); % half times r*Cos(45deg) magnified, divided by pixel size
                    size_central_insert = 5*magnification/pixel_size_1;
                    size_tungsten = 1*magnification/(pixel_size_1);

                    ROIsCal = CalcROIs(centresCalcium, sizes);
                    ROIsIod = CalcROIs(centresIodine, sizes);
                    ROIsMixed = CalcROIs(centresMixed, sizes);
                    ROIsTungsten = CalcROIs(centresTungsten, size_tungsten);
                    ROIcentral_insert = CalcROIs([centre_p,centre_p], size_central_insert);
                    ROIcentral_MTF = CalcROIs([centre_p,centre_p], size_central_insert*2);

                    waterROI_offset = 1.5*size_central_insert;
                    centresWater = [centre_p+waterROI_offset,centre_p+waterROI_offset;
                        centre_p+waterROI_offset,centre_p-waterROI_offset;
                        centre_p-waterROI_offset,centre_p+waterROI_offset;
                        centre_p-waterROI_offset,centre_p-waterROI_offset;];
                    ROIsWater = CalcROIs(centresWater, size_central_insert);
                end

                ROI_CNRandSNR = [int16(ROIcentral_insert); int16(ROIsWater(1,:))];
                ROIs_NPS = int16(ROIsWater);
                ROIs_Iodine = int16(ROIsIod);
                ROIs_Calcium = int16(ROIsCal);

                for DisplayROIs = 1:1

                    %Display ROIs
                    figure(2);
                    imshow(ImageData,[]);
                    title('Image with ROIs');

                    PlotROIs(ROI_CNRandSNR, 2, 'rectangle','b', 2)
                    %PlotROIs(ROI_MTF, 2, 'circle', 'r', 2)
                    PlotROIs(ROIcentral_insert, 2, 'circle', 'r', 2)
                    PlotROIs(ROIcentral_MTF, 2, 'circle', 'r', 2)
                    PlotROIs(ROIs_NPS, 2, 'rectangle', 'g', 2)
                    PlotROIs(ROIs_Iodine, 2, 'rectangle', 'm', 2)
                    PlotROIs(ROIs_Calcium, 2, 'rectangle', 'w', 2)
                    %PlotROIs(ROI_Mixtures, 2, 'rectangle', 'c', 2)

                end



                %calculate MTF

                plotMTF = 0; %Change to 1 if you want to see graphics associated with the progressive calculation of MTF
                %MTF(Flux,Pitch,Alg,Approach,:) = MTF_Calculator(ImageData, ROI_MTF, plotMTF, 3);



                %calculate CNR and SNR

                CNR(iFlux,iPitch,Alg,recon) = CalcCNR(ImageData, ROI_CNRandSNR)

                SNR(iFlux,iPitch,Alg,recon) = CalcSNR(ImageData, ROI_CNRandSNR);



                %Calculate NPS

                plotNPS = 1; %Change to 1 if you want to see 2D and 1D plots of the calculated NPS
                 newNPS = Calculate_NPS(ImageData, ROIs_NPS, plotNPS, 4);
                    NPS(iFlux,iPitch,Alg,recon,1:size(newNPS,1),1:size(newNPS,2)) = newNPS;


                %Calculate linearity of I conc

                plotILin = 1; %Change to 1 if you want to see linearity plot of I concentration vs Measured signal
                [IodineSignals(iFlux,iPitch,Alg,recon,:), IodineLinearity(iFlux,iPitch,Alg,recon)] = LinearityCalculator(ImageData,ROIs_Iodine, IodineConcentrations_Actual, plotILin, 5);



                %Calculate linearity of Ca conc (can be same function as previous)

                plotCaLin = 1; %Change to 1 if you want to see linearity plot of Ca concentration vs Measured signal
                [CalciumSignals(iFlux,iPitch,Alg,recon,:), CalciumLinearity(iFlux,iPitch,Alg,recon)] = LinearityCalculator(ImageData,ROIs_Calcium, CalciumConcentrations_Actual, plotCaLin, 6);



                %Predict I and Ca concentrations in unknown ROIs...

                %...And compare it to the expected/known values


               pause(1);
            end
        end
    end
end

SAVEOVER = 0;
if SAVEOVER == 1
save('IQA_statsArrays_t1','CNR','SNR','NPS',"CalciumLinearity","CalciumSignals","IodineLinearity","IodineSignals",'-V7.3');
end

%%
ANALYSEOLD = 1;
if ANALYSEOLD == 1

end

%%
clc;

act_idxs = 3:6;     % approaches 3–6 are ACTs
add_idxs = 7:12;    % approaches 7–12 are additive CSCAs
sub_idxs = 13:18;   % approaches 13–18 are subtractive CSCAs

CNR_std = CNR(:,:,1);
SNR_std = SNR(:,:,1);

[CNR_act, CNR_max_act, CNR_maxFluxIdx_act, CNR_maxPitchIdx_act, CNR_maxAlgIdx_act] = maxFinder(CNR,act_idxs);
[CNR_add, CNR_max_add, CNR_maxFluxIdx_add, CNR_maxPitchIdx_add, CNR_maxAlgIdx_add] = maxFinder(CNR,add_idxs);
[CNR_sub, CNR_max_sub, CNR_maxFluxIdx_sub, CNR_maxPitchIdx_sub, CNR_maxAlgIdx_sub] = maxFinder(CNR,sub_idxs);
[SNR_act, SNR_max_act, SNR_maxFluxIdx_act, SNR_maxPitchIdx_act, SNR_maxAlgIdx_act] = maxFinder(SNR,act_idxs);
[SNR_add, SNR_max_add, SNR_maxFluxIdx_add, SNR_maxPitchIdx_add, SNR_maxAlgIdx_add] = maxFinder(SNR,add_idxs);
[SNR_sub, SNR_max_sub, SNR_maxFluxIdx_sub, SNR_maxPitchIdx_sub, SNR_maxAlgIdx_sub] = maxFinder(SNR,sub_idxs);

%%


function [slices, mValue, mFlux, mPitch, mAlg] = maxFinder(array,indices);
% Get max for each class
for iFlux = 1:3  
    arraySlice = array(iFlux,:,indices);
    slices(iFlux,:,:) = arraySlice;
    mValue(iFlux) = max(arraySlice(:));
    best_idx_act(iFlux) = find(array == mValue(iFlux),1); %%This returns first value only. Upon testing, performance pairs are 1 - 5 (STD and PCS) and 3 - 6 (DLR and SR). These will only differ at higher thresholds.
    [mFlux(iFlux),mPitch(iFlux),mAlg(iFlux)] = ind2sub(size(array),best_idx_act(iFlux));

end
end



%%
%{


clear
clc
close all
load('LocationsForCERN2024Images','RECT_CNR_B','RECT_CNR_C','RECT_IConc','RECT_CaConc','RECT_MixConc')

boxsize = 16;

%READ IN IMAGE WITH 0 THRESHOLD LABEL
%imagepath = 'C:\Users\oliep\Downloads\results sharing (3)\LF_200\threshold';
%imagepath = 'C:\Users\oliep\OneDrive\Work\dataForCern2024\100um\threshold_0\';
%imagename = '\filtered_backprojection_Combined_sinograms_ADD2x2Dy_threshold_';
imagename = 'backprojection_Combined_sinograms_Add2x2Dy_threshold_0.png';
%filelocation = strcat(imagepath,num2str(0),imagename,num2str(0),'.png')
%filelocation = imagename;

imagename = 'C:\Users\oliep\backprojection_Combined_sinograms_Add2x2Dy_threshold_0.png'
filelocation = imagename;

uiopen(filelocation,true)
pause(13)
restarting = 1%%

image = cdata(:,:,1);
%%
% Display the image
figure(1);
imshow(image);
imcontrast();
title('Select an area by drawing a square');
pause(7)
restarting
%%

%get contrast area
%RECT_CNR_B = getrect;
xCNR_B = RECT_CNR_B(1);
yCNR_B = RECT_CNR_B(2);
widthCNR = boxsize;%RECT_CNR_B(3);
heightCNR = boxsize;%RECT_CNR_B(4);
rectangle('Position', [xCNR_B, yCNR_B, widthCNR, heightCNR], 'EdgeColor', 'c', 'LineWidth', 2);
ROIB = image(yCNR_B:(yCNR_B + heightCNR), xCNR_B:(xCNR_B + widthCNR));
ROIBackground(1) = mean2(ROIB);

%get background area
%RECT_CNR_C = getrect;
xCNR_C = RECT_CNR_C(1);
yCNR_C = RECT_CNR_C(2);
rectangle('Position', [xCNR_C, yCNR_C, widthCNR, heightCNR], 'EdgeColor', 'm', 'LineWidth', 2);
ROIC = image(yCNR_C:(yCNR_C + heightCNR), xCNR_C:(xCNR_C + widthCNR));
ROIContrast(1) = mean2(ROIC);

%Calculate CNR
STDBackground(1) = std(double(ROIB),0,"all");
CNR(1) = (ROIContrast(1) - ROIBackground(1))/STDBackground(1)





for iijjkk = 1:6
    %RECT_CaConc(:,iijjkk) = getrect;
    x_I = RECT_IConc(1,iijjkk);
    y_I = RECT_IConc(2,iijjkk);
    widthI = boxsize%RECT_IConc(3,1);
    heightI = boxsize%RECT_IConc(4,1);
    switch iijjkk
        case 1
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'r', 'LineWidth', 2);
        case 2
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'c', 'LineWidth', 2);
        case 3
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'y', 'LineWidth', 2);
        case 4
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'g', 'LineWidth', 2);
        case 5
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'b', 'LineWidth', 2);
        case 6
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'k', 'LineWidth', 2);
    end
    ROII = image(y_I-heightI/2:(y_I + heightI/2-1), x_I-widthI/2:(x_I + widthI/2-1));
    ConcI(iijjkk,1) = mean2(ROII);

end





for iijjkk = 1:6
    %RECT_CaConc(:,iijjkk) = getrect;
    x_I = RECT_CaConc(1,iijjkk);
    y_I = RECT_CaConc(2,iijjkk);
    widthI = boxsize%RECT_IConc(3,1);
    heightI = boxsize%RECT_IConc(4,1);
    switch iijjkk
        case 1
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'r', 'LineWidth', 2);
        case 2
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'c', 'LineWidth', 2);
        case 3
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'y', 'LineWidth', 2);
        case 4
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'g', 'LineWidth', 2);
        case 5
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'b', 'LineWidth', 2);
        case 6
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'k', 'LineWidth', 2);
    end
    ROII = image(y_I-heightI/2:(y_I + heightI/2-1), x_I-widthI/2:(x_I + widthI/2-1));
    ConcCa(iijjkk,1) = mean2(ROII);

end






for iijjkk = 1:6
    %RECT_CaConc(:,iijjkk) = getrect;
    x_I = RECT_MixConc(1,iijjkk);
    y_I = RECT_MixConc(2,iijjkk);
    widthI = boxsize%RECT_IConc(3,1);
    heightI = boxsize%RECT_IConc(4,1);
    switch iijjkk
        case 1
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'r', 'LineWidth', 2);
        case 2
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'c', 'LineWidth', 2);
        case 3
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'y', 'LineWidth', 2);
        case 4
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'g', 'LineWidth', 2);
        case 5
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'b', 'LineWidth', 2);
        case 6
            rectangle('Position', [x_I-widthI/2, y_I-heightI/2,  widthI, heightI], 'EdgeColor', 'k', 'LineWidth', 2);
    end
    ROII = image(y_I-heightI/2:(y_I + heightI/2-1), x_I-widthI/2:(x_I + widthI/2-1));
    ConcMix(iijjkk,1) = mean2(ROII);

end

labell = 11;
[MI(1), CI(1), RI(1)] = ScatterPlotter(ConcI(:,1), labell);
movegui('west')
labell = 21;
[MCa(1), CCa(1), RCa(1)] = ScatterPlotter(ConcCa(:,1), labell);
movegui('east')



%%

pause(3)
for ijk = 1:3
    filelocation = strcat(imagepath,num2str(ijk),imagename,num2str(ijk),'.png')
    clear cdata
    uiopen(filelocation,1)
    pause(3)
    image = cdata(:,:,1);
    figure(ijk+1)
    imshow(image)
    %get contrast area
    rectangle('Position', [xCNR_C, yCNR_C, widthCNR, heightCNR], 'EdgeColor', 'c', 'LineWidth', 2);
    ROIC = image(yCNR_C:(yCNR_C + heightCNR), xCNR_C:(xCNR_C + widthCNR));
    ROIContrast(ijk+1) = mean2(ROIC);

    %get background area
    rectangle('Position', [xCNR_B, yCNR_B, widthCNR, heightCNR], 'EdgeColor', 'm', 'LineWidth', 2);
    ROIB = image(yCNR_B:(yCNR_B + heightCNR), xCNR_B:(xCNR_B + widthCNR));
    ROIBackground(ijk+1) = mean2(ROIB);

    %Calculate CNR
    STDBackground(ijk+1) = std(double(ROIB),0,"all");
    CNR(ijk+1) = (ROIContrast(ijk+1) - ROIBackground(ijk+1))/STDBackground(ijk+1)



    for iijjkk = 1:6
        x_I = RECT_IConc(1,iijjkk);
        y_I = RECT_IConc(2,iijjkk);
        %widthI = RECT_IConc(3,1);
        %heightI = RECT_IConc(4,1);
        switch iijjkk
            case 1
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'r', 'LineWidth', 2);
            case 2
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'c', 'LineWidth', 2);
            case 3
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'y', 'LineWidth', 2);
            case 4
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'g', 'LineWidth', 2);
            case 5
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'b', 'LineWidth', 2);
            case 6
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'k', 'LineWidth', 2);
        end
        ROII = image(y_I-heightI/2:(y_I + heightI/2-1), x_I-widthI/2:(x_I + widthI/2-1));
        ConcI(iijjkk,ijk+1) = mean2(ROII);
    end



    for iijjkk = 1:6
        x_I = RECT_CaConc(1,iijjkk);
        y_I = RECT_CaConc(2,iijjkk);
        %widthI = RECT_IConc(3,1);
        %heightI = RECT_IConc(4,1);
        switch iijjkk
            case 1
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'r', 'LineWidth', 2);
            case 2
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'c', 'LineWidth', 2);
            case 3
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'y', 'LineWidth', 2);
            case 4
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'g', 'LineWidth', 2);
            case 5
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'b', 'LineWidth', 2);
            case 6
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'k', 'LineWidth', 2);
        end
        ROII = image(y_I-heightI/2:(y_I + heightI/2-1), x_I-widthI/2:(x_I + widthI/2-1));
        ConcCa(iijjkk,ijk+1) = mean2(ROII);
    end

    for iijjkk = 1:6
        x_I = RECT_MixConc(1,iijjkk);
        y_I = RECT_MixConc(2,iijjkk);
        %widthI = RECT_IConc(3,1);
        %heightI = RECT_IConc(4,1);
        switch iijjkk
            case 1
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'r', 'LineWidth', 2);
            case 2
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'c', 'LineWidth', 2);
            case 3
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'y', 'LineWidth', 2);
            case 4
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'g', 'LineWidth', 2);
            case 5
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'b', 'LineWidth', 2);
            case 6
                rectangle('Position', [x_I-widthI/2, y_I-heightI/2, widthI, heightI], 'EdgeColor', 'k', 'LineWidth', 2);
        end
        ROII = image(y_I-heightI/2:(y_I + heightI/2-1), x_I-widthI/2:(x_I + widthI/2-1));
        ConcMix(iijjkk,ijk+1) = mean2(ROII);

    end




    labell = 11 + ijk;
    [MI(ijk+1), CI(ijk+1), RI(ijk+1)] =  ScatterPlotter(ConcI(:,ijk+1), labell);
    movegui('west')

    labell = 21 + ijk;
    [MCa(ijk+1), CCa(ijk+1), RCa(ijk+1)] =  ScatterPlotter(ConcCa(:,ijk+1), labell);
    movegui('east')
end



labell = 14;
Conc_unweighted = sum(ConcI(:,:),2)/4
[foo, bar, RI(5)] = ScatterPlotter(Conc_unweighted,labell+1);
Conc_weighted = (ConcI(:,1)*CNR(1) + ConcI(:,2)*CNR(2) + ConcI(:,3)*CNR(3) + ConcI(:,4)*CNR(4)) / sum(CNR)
ScatterPlotter(Conc_weighted,labell+2)

labell = 24;
Conc_unweighted = sum(ConcCa(:,:),2)/4
[foo2, bar2, RI(6)] = ScatterPlotter(Conc_unweighted,labell+1);
Conc_weighted = (ConcCa(:,1)*CNR(1) + ConcCa(:,2)*CNR(2) + ConcCa(:,3)*CNR(3) + ConcCa(:,4)*CNR(4)) / sum(CNR)
ScatterPlotter(Conc_weighted,labell+2)

figure(11)
movegui('northwest')
figure(12)
movegui('southwest')
figure(13)
movegui('north')
figure(14)
movegui('south')
figure(15)
movegui('northeast')
figure(16)
movegui('southeast')

figure(21)
movegui('northwest')
figure(22)
movegui('southwest')
figure(23)
movegui('north')
figure(24)
movegui('south')
figure(25)
movegui('northeast')
figure(26)
movegui('southeast')
%% Calculate I concentrations

%coefficients(1), coefficients(2)
%MI(1:4) = [0.04489, 0.03156, 0.02347, 0.02677];
%CI(1:4) = [89.03674, 96.65997, 104.46972, 97.99927];


for kji = 1:4
    Estimate(:,kji) = (ConcMix(:,kji) - CI(kji))/MI(kji);
end
Estimate


MWAv = foo;
CWAv = bar;

Estimate2(:,1) = (sum(ConcMix(:,:),2)/4 - CWAv/MWAv);

Estimate2

%{
figure(111)
hold off
LayeredPlotter([100, 80, 60, 40, 20, 0],Estimate(:,1),111,'r')
hold on
LayeredPlotter([100, 80, 60, 40, 20, 0],Estimate(:,2),111,'c')
LayeredPlotter([100, 80, 60, 40, 20, 0],Estimate(:,3),111,'g')
LayeredPlotter([100, 80, 60, 40, 20, 0],Estimate(:,4),111,'b')
hold off
%}

%%
%{
for HF250 = 1:1
%Calculate I concentrations
MI(1:4) = [0.06, 0.11, 0.17, 0.19];
CI(1:4) = [130.14, 135.56, 122.16, 115.67];

%Calculate Ca concentrations

MCa(1:4) = [0.00, 0.01, 0.02, 0.02];
CCa(1:4) = [130.21, 135.64, 122.23, 115.70];
end
%}

%{
for HF100_NOSTD = 1:1
%Calculate I concentrations
MI(1:4) = [0.03932, 0.03194, 0.02253, 0.03234];
CI(1:4) = [102.37828, 96.71128, 103.72646, 105.31217];

%Calculate Ca concentrations

MCa(1:4) = [0.00584, 0.00309, 0.00360, 0.00365];
CCa(1:4) = [102.31189, 96.76870, 103.66912, 105.27293];
end
%}
%{
for HF100_NOSR = 1:1
%Calculate I concentrations
MI(1:4) = [0.04489, 0.03156, 0.02347, 0.02677];
CI(1:4) = [89.03674, 96.65997, 104.46972, 97.99927];

%Calculate Ca concentrations

MCa(1:4) = [0.00307, 0.00241, 0.00303, 0.00563];
CCa(1:4) = [89.04688, 96.71036, 104.52077, 97.94261];
end
%}

x = [10, 28, 41, 55];

SpectrumA = 100*MI(:) + CI(:); %ConcI(1,:);
SpectrumB = 100*MCa(:) + CCa(:); %ConcCa(1,:);

for jki = 1:6
    clear fitParams
    mixedData = ConcMix(jki,:);

    modelFunction = @(theta,x) min(max(theta(1), 0), 1) * SpectrumA' + min(max(theta(2), 0), 1) * SpectrumB';

    initialGuess = [0.01,0.99];

    options = optimset('Display', 'iter');

    fitParams = lsqcurvefit(modelFunction, initialGuess, x, mixedData, [], [], options);

    thetaA(jki) = min(max(fitParams(1), 0 ),1);
    thetaB(jki) = min(max(fitParams(2),0),1);

    figure(3748)
    hold off
    plot(x,thetaA(jki)*SpectrumA+thetaB(jki)*SpectrumB, '--or');
    hold on
    plot(x,mixedData,'-.ok')
    plot(x,SpectrumA,'-og')
    plot(x,SpectrumB,'-ob')
    title('Estimation of I and Ca components')
    xlabel('Threshold number')
    ylabel('Pixel value')
    legend('Estimated combination spectrum','Measured spectrum','Iodine spectrum','Calcium spectrum', 'location','southwest')
    pause(1)
end

thetaA
thetaB





%%
function [COEFFM, COEFFC, R2] = ScatterPlotter(CONC,labell)
ConcAscending(1:6) = 0;
for ijk = 1:6
    ConcAscending(ijk) = CONC(7-ijk);
end
figure(labell);
x = [0.5 1 5 10 50 100];
y = ConcAscending;

scatter(x, y, 'b');  % Scatter plot
hold on;

coefficients = polyfit(x, y, 1);  % Fit a first-degree polynomial (line)

% Compute the fitted values for the line
numFitPoints = 1000;  % Enough points to make the plot look continuous
xFit = linspace(min(x), max(x), numFitPoints);
yFit = polyval(coefficients, xFit);

% Plot the line of best fit
plot(xFit, yFit, 'r-', 'LineWidth', 2);
xlabel('Iodine concentration (mg/ml)');
ylabel('Average pixel value');
title('Iodine quantification');
COEFFM = coefficients(1);
COEFFC = coefficients(2);
% Display the equation of the line
equation = sprintf('y = %.5fx + %.5f', coefficients(1), coefficients(2));
text(max(x), max(y), equation, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Calculate R² (coefficient of determination)
yPredicted = polyval(coefficients, x);
SSres = sum((y - yPredicted).^2);
SStot = sum((y - mean(y)).^2);
rSquared = 1 - SSres / SStot;
rSquaredText = sprintf('R² = %.4f', rSquared);
text(max(x), max(y) - 0.1 * (max(y) - min(y)), rSquaredText, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

hold off;
R2 = rSquared;
end

function LayeredPlotter(x,y,fignum,colour)
figure(fignum)
scatter(x,y,colour)
coefficients = polyfit(x, y, 1);
% Compute the fitted values for the line
numFitPoints = 1000;  % Enough points to make the plot look continuous
xFit = linspace(min(x), max(x), numFitPoints);
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, colour, 'LineWidth', 2);

xlabel('X values');
ylabel('Y values');
title('Scatter Plot with Line of Best Fit');

% Display the equation of the line
equation = sprintf('y = %.2fx + %.2f', coefficients(1), coefficients(2))
%text(max(x), max(y), equation, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Calculate R² (coefficient of determination)
yPredicted = polyval(coefficients, x);
SSres = sum((y - yPredicted).^2);
SStot = sum((y - mean(y)).^2);
rSquared = 1 - SSres / SStot;
rSquaredText = sprintf('R² = %.4f', rSquared)
%text(max(x), max(y) - 0.1 * (max(y) - min(y)), rSquaredText, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

end

%}



%% FUNCTIONS

function [ImageData] = LoadAndProcessImage(ImageFileName)
%loads in the given image
rawImage = imread(ImageFileName);

% Check if the image is color and convert to grayscale if necessary
if size(rawImage, 3) == 3
    % Convert to grayscale
    grayscale_image = rgb2gray(rawImage);
else
    % Directly assign grayscale image
    grayscale_image = ct_image_uint;
end

% Convert to double precision for best accuracy
ImageData = double(grayscale_image);
end



function [ROIs] = CalcROIs(centres,sizes)
%Produces a 4 column array containing x_lower, x_upper, y_lower and y_upper
%for each of the listed centres. Lower (and Upper) given by centre minus
%(and plus) half of the relevant ROI size value.

halfLengths = round(sizes/2);

ROIs = [centres(:,1) - halfLengths, centres(:,1) + halfLengths, ... % New columns 1 and 2
    centres(:,2) - halfLengths, centres(:,2) + halfLengths];    % New columns 3 and 4

end



function [ROIs] = TranslateROIs(refROI,refPitch,currentPitch)

refCentres = [ ( refROI(:,1)+refROI(:,2) )/2, ( refROI(:,3)+refROI(:,4) )/2 ];
refSizes = refROI(:,2) - refROI(:,1);

ScalingRatio = currentPitch / refPitch;

ScaledSizes = (refSizes * ScalingRatio);
ScaledCentres = round(refCentres * ScalingRatio);

ROIs = CalcROIs(ScaledCentres,ScaledSizes);

end



function PlotROIs(ROIs, ImageHandle, shape, c,w)

figure(ImageHandle)
hold on;

if shape == "rectangle"
    % Draw a rectangle around the phantom ROIs
    for ij = 1:size(ROIs,1)
        fill([ROIs(ij,1), ROIs(ij,2), ROIs(ij,2), ROIs(ij,1)], [ROIs(ij,4), ROIs(ij,4), ROIs(ij,3), ROIs(ij,3)], 'r', 'FaceColor', 'none', 'EdgeColor', c, 'LineWidth', w);
    end
elseif shape == "circle"
    % Draw a circle around the phantom ROIs
    for ij = 1:size(ROIs,1)
        x_centre = ( ROIs(ij,2) + ROIs(ij,1) ) / 2;
        y_centre = ( ROIs(ij,4) + ROIs(ij,3) ) / 2;
        circle_centre = [x_centre, y_centre];
        circle_radius = (ROIs(ij,2) - ROIs(ij,1) ) / 2;
        % Plot the circle
        viscircles(circle_centre, circle_radius, 'EdgeColor', c, 'LineWidth', w);
    end
end

hold off

end



function MTF = MTF_Calculator(imageData, ROI, PLOT, ImageHandle)

%Identify circle
inputCircleCentre = center
inputCircleRadius = radius
inputOuterCircleMultiplier = outer_radius_multiplier

%Identify regions for data normalisation
%ROI needed for innermaterial square = innerMaterial_ROIData
%ROI needed for outermaterial square = outerMaterial_ROIData

%normalise data based on inner and outer material data
avg_inner = mean2(innerMaterial_ROIData);
avg_outer = mean2(outerMaterial_ROIData);


% Apply the normalization formula
normalisedImageData = (imageData - avg_outer) / (avg_inner - avg_outer);


%Calculate ESF
max_radius = radius*outer_radius_multiplier;
viscircles(center, max_radius,'EdgeColor','b');


%Make radial lines to outside of ring
angularStepSize = 1;
theta_degrees = 0:angularStepSize:360-angularStepSize; % Angles in degrees
theta_radians = deg2rad(theta_degrees); % Convert to radians

% Convert polar coordinates to Cartesian coordinates
SamplingFrequency = 1
intensityValues(1:360/angularStepSize,1:uint16(max_radius*SamplingFrequency)) = 0;
if PLOT == 1
    figure(ImageHandle*100+4)
    hold off
    figure(ImageHandle*100+5)
    hold off
    movegui("northwest")
    figure(ImageHandle*100+6)
    hold off
    movegui("southwest")
end


for rho = 0:1/SamplingFrequency:max_radius
    rho
    [x, y] = pol2cart(theta_radians, rho);

    %work out radius in pol coordinates
    figure(ImageHandle*100+4)
    if PLOT == 1
        imshow(normalisedImageData)
        hold on
        viscircles(center, radius,'EdgeColor','r');
        viscircles(center, max_radius,'EdgeColor','b');
        viscircles(center, rho,'EdgeColor','g');
    end
    movegui("center")

    % Sample intensity values
    intensityValues(1:360/angularStepSize,uint16(rho*SamplingFrequency)+1) = interp2(normalisedImageData, center(1) +x, center(2)+y, 'linear');
    if PLOT == 1
        figure(ImageHandle*100+5)
        plot(intensityValues(:,uint16(rho*SamplingFrequency)))
        hold on
        title(rho*SamplingFrequency)
        figure(ImageHandle*100+6)
        polarplot(theta_radians, intensityValues);
        title('Intensity Profile Along Radial Lines');
        pause(0.01)
    end

end


if PLOT == 1
    figure(ImageHandle*100+5)
    movegui("northeast")
    for ii = 2:360/angularStepSize
        %hold off
        %for jj = 1:ii-1
        plot(1:size(intensityValues,2),intensityValues(ii-1,:),'k')
        hold on
        %end
        plot(1:size(intensityValues,2),intensityValues(ii,:),'r')
        %hold off
        title(ii)
        pause(0.05)
    end

    figure(ImageHandle*100+6)
    movegui("southeast")
    plot(mean(intensityValues(:,:),1))
end

%%
%Differentiate ESF to get LSF
% this is the "finite difference" derivative. Note it is  one element shorter than y and x
clear x
x = [1:length(mean(intensityValues(:,:),1))];
yd = (diff(mean(intensityValues(:,:),1))./diff(x));
LSF = yd;
% this is to assign yd an abscissa midway between two subsequent x
xd = (x(2:end)+x(1:(end-1)))/2;
% this should be a rough plot of derivative
if PLOT == 1
    figure(7)
    hold off
    plot(xd,yd)
end



L = length(LSF); % Length of the LSF

%%If Hann filtering
hann_window = hann(L); % Generate the Hann window
LSF_hann = LSF .* hann_window'; % Apply the window
LSF = LSF_hann;
%%If Hann filtering

if PLOT == 1
    figure(7)
    hold on
    plot(xd,LSF_hann)
    legend('raw','hann-filtered')
    hold off
end

P2 = abs(fft(LSF / L)); % Compute the FFT
P1 = P2(1:floor(L/2)+1); % Take only the positive frequencies
%P1(2:end-1) = 2 * P1(2:end-1); % Double the amplitudes (except DC and Nyquist)

nP1 = P1 - min(P1); % Shift to non-negative values
nP1 = nP1 ./ max(nP1); % Normalize to [0, 1]

Fs = SamplingFrequency; % Sampling frequency (pixels/lp)
f = Fs * (0:(L/2)) / L; % Normalized spatial frequency

if PLOT == 1
    figure(8)
    hold off
    plot(f(1:end), nP1(1:end), 'o-','linewidth', 3);
    title('Radial MTF');
    xlabel('Normalized Spatial Frequency');
    ylabel('MTF');
end





end



function [CNR] = CalcCNR(imageData, ROIs)

%%This funciton assumes that the first ROI value will be the contrast agent and
%%the second will be the bakcground

ContrastData = imageData( ROIs(1,1):ROIs(1,2) , ROIs(1,3):ROIs(1,4) );
BackgroundData = imageData( ROIs(2,1):ROIs(2,2) , ROIs(2,3):ROIs(2,4) );

%Calculate CNR
ContrastMean = mean2(ContrastData);
BackgroundMean = mean2(BackgroundData);
BackgroundStd = std(BackgroundData,0,"all");

CNR = (ContrastMean - BackgroundMean)/BackgroundStd;

end



function [SNR] = CalcSNR(imageData, ROIs)

%%This funciton assumes that the first ROI value will be the Signal and
%%the second will be the background. ALSO NOTE that this current
%%implementation assumes that the CNR and SNR use the same ROIs, i.e.
%%Contrast agent = Signal, background = homogenous phantom region without
%%contrast agent in it.

SignalData = imageData( ROIs(1,1):ROIs(1,2) , ROIs(1,3):ROIs(1,4) );
BackgroundData = imageData( ROIs(2,1):ROIs(2,2) , ROIs(2,3):ROIs(2,4) );

%Calculate SNR
SignalMean = mean2(SignalData);
%SignalStd = std(SignalData,0,"all");
BackgroundStd = std(BackgroundData,0,"all");

%SNR = SignalMean / SignalStd;
SNR = SignalMean / BackgroundStd;

end



function [NPS] = Calculate_NPS(imageData, ROIs, PLOT, ImageHandle)
% Function to calculate NPS for a single ROI in a given image.

% Extract the current ROI
currentROIData = imageData( ROIs(1,1):ROIs(1,2) , ROIs(1,3):ROIs(1,4) );

% Subtract mean from ROI to isolate noise
currentROINoiseEst = currentROIData - mean(currentROIData(:));

% Perform 2D Fourier Transform
fft_result = fft2(currentROINoiseEst);
fft_shifted = fftshift(fft_result);

% Calculate Power Spectrum
noisePowerSpectrum = abs(fft_shifted).^2;

% Normalize by ROI size
[rows, cols] = size(ROIs);
NPS = noisePowerSpectrum / (rows * cols);

% Display the Noise Power Spectrum
if PLOT == 1
figure(ImageHandle);
imagesc(log(NPS + 1)); % Use log scale for better visualization
colormap('jet');
colorbar;
title('Noise Power Spectrum (NPS)');

    % Optionally extract the 1D NPS along specific directions, e.g., radial averaging
    [fx, fy] = meshgrid(linspace(0, 50, cols), linspace(0, 50, rows)); % Frequency axes
    radial_distances = sqrt(fx.^2 + fy.^2);
    unique_distances = unique(radial_distances(:));
    nps_1d = zeros(size(unique_distances));

    for i = 1:length(unique_distances)
        mask = abs(radial_distances - unique_distances(i)) < 1e-3; % Small threshold for matching
        nps_1d(i) = mean(NPS(mask));
    end

    % Plot the 1D Noise Power Spectrum
    figure(ImageHandle*100);
    plot(unique_distances, nps_1d, '-o', 'LineWidth', 2);
    title('1D Noise Power Spectrum');
    xlabel('Spatial Frequency');
    ylabel('Power');
    grid on;
end
end



function [CConc_Sig, Linearity] = LinearityCalculator(imageData, ROIs, CConc_Act, PLOT, ImageHandle)
CConc_Sig = zeros(size(ROIs, 1), 1);
for ij = 1:size(ROIs,1)
    x1 = ROIs(ij,1);
    x2 = ROIs(ij,2);
    y1 = ROIs(ij,3);
    y2 = ROIs(ij,4);
    CConc_Sig(ij) = mean2(imageData(y1:y2,x1:x2));
end
%figure(11)
%imshow(imageData,[])

CConc_Act; % Material concentration (column 1)
CConc_Act = CConc_Act';
CConc_Sig; % Measured signal (column 2)

coefficients = polyfit(CConc_Act, CConc_Sig, 1); % Linear fit (1st degree polynomial)

if PLOT == 1
    % Compute the fitted values for the line
    numFitPoints = (size(CConc_Act) + 2) * 10;  % Enough points to make the plot look continuous
    maxval = max(CConc_Act);
    minval = min(CConc_Act);
    numGivenPoints = size(CConc_Act,1);
    deltaval = (maxval - minval)/numGivenPoints;
    numFitPoints = (numGivenPoints + 2) * 10;
    xFit = linspace(minval - deltaval,  maxval + deltaval, numFitPoints); %This line determines
    yFit = polyval(coefficients, xFit);


    figure(ImageHandle)
    plot(CConc_Act, CConc_Sig, 'o', 'DisplayName', 'Data Points'); % Scatter plot of data
    hold on;
    plot(xFit, yFit, '-r', 'DisplayName', 'Fitted Line'); % Fitted line
    xlabel('Material Concentration');
    ylabel('Measured Signal');
    title('Contrast concentration vs Signal intensity');
    Coeff_M = coefficients(1);
    Coeff_C = coefficients(2);
    legend show;
    grid on;
    hold off
    % Display the equation of the line
    equation = sprintf('y = %.5fx + %.5f', Coeff_M, Coeff_C);
    text(max(CConc_Act), max(CConc_Sig), equation, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end

% Calculate R² (coefficient of determination)
yPredicted = polyval(coefficients, CConc_Act);
y_mean = mean(CConc_Sig); % Mean of the measured signal
SStot = sum((CConc_Sig - y_mean).^2); % Total sum of squares
SSres = sum((CConc_Sig - yPredicted).^2);  % Residual sum of squares
R2 = 1 - (SSres / SStot);     % Coefficient of determination

%disp(['R-squared: ', num2str(R2)]);

Linearity = R2;

end






%{

%%NPS V 1-
    % Extract the current ROI
    roiCoords = rois{i};
    x = roiCoords(1);
    y = roiCoords(2);
    width = roiCoords(3);
    height = roiCoords(4);
    roi = image(y:y+height-1, x:x+width-1);

    % Subtract mean from ROI to isolate noise
    roi_mean_subtracted = roi - mean(roi(:));

    % Perform 2D Fourier Transform
    fft_result = fft2(roi_mean_subtracted);
    fft_shifted = fftshift(fft_result);

    % Calculate Power Spectrum
    power_spectrum = abs(fft_shifted).^2;

    % Normalize by ROI size
    nps = power_spectrum / (width * height);

    % Store individual NPS
    npsList{i} = nps;

    % Accumulate for average calculation
    totalNPS = totalNPS + nps;


%%    NPS - V2
%% Compute the average NPS across all ROIs
avgNPS = totalNPS / length(rois);

% Display the average NPS
figure;
imagesc(log(avgNPS + 1)); % Log scale for better visualization
colormap('jet');
colorbar;
title('Average Noise Power Spectrum');
end



% Ensure the ROI is in double format for calculations
ROI = double(ROI);

% Subtract the mean intensity from the ROI to isolate the noise
roi_mean_subtracted = ROI - mean(ROI(:));

% Perform a 2D Fourier Transform on the ROI
fft_result = fft2(roi_mean_subtracted);
fft_shifted = fftshift(fft_result); % Shift the zero frequency component to
the centre

% Calculate the Power Spectrum
power_spectrum = abs(fft_shifted).^2;

% Normalize by the area (number of pixels) in the ROI
[rows, cols] = size(ROI);
nps = power_spectrum / (rows * cols);

% Display the Noise Power Spectrum
figure;
imagesc(log(nps + 1)); % Use log scale for better visualization
colormap('jet');
colorbar;
title('Noise Power Spectrum (NPS)');

% Optionally extract the 1D NPS along specific directions, e.g., radial averaging
[fx, fy] = meshgrid(linspace(-0.5, 0.5, cols), linspace(-0.5, 0.5, rows)); % Frequency axes
radial_distances = sqrt(fx.^2 + fy.^2);
unique_distances = unique(radial_distances(:));
nps_1d = zeros(size(unique_distances));

for i = 1:length(unique_distances)
    mask = abs(radial_distances - unique_distances(i)) < 1e-3; % Small threshold for matching
    nps_1d(i) = mean(nps(mask));
end

% Plot the 1D Noise Power Spectrum
figure;
plot(unique_distances, nps_1d, 'LineWidth', 2);
title('1D Noise Power Spectrum');
xlabel('Spatial Frequency');
ylabel('Power');
grid on;

%}