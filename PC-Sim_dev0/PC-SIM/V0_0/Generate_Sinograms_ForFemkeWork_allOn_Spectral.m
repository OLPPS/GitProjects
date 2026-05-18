%make a sinogram
clear
clc
BASEFILENAME = 'FullPhantom_Smaller_RealisticContrast';
IDSTRING = '10e7';
XMax = 292;
NumberOfProjections_Total = 2000;
NumberOfProjections_PerFile = 40;
NumOfThresh = 4;
NumberOfRuns = 50
NumOfApproaches = 18;
sinograms = zeros(XMax,NumberOfProjections_Total,NumOfThresh,NumOfApproaches);

for RUN_number = 0:NumberOfRuns - 1
    for projection =  1:NumberOfProjections_PerFile
    inputfilename = strcat(BASEFILENAME,'_',IDSTRING,'_','Run_',num2str(RUN_number),'_projection_',num2str(projection),'_',num2str(NumOfThresh),'Thresh_projectionData.mat')
    load(inputfilename);
    sinogramColumn = RUN_number*NumberOfProjections_PerFile + projection

    sinograms(:,sinogramColumn,:,1) = NoCSCA_STD_Counters(:,1,:);
    
    sinograms(:,sinogramColumn,:,2) = NoCSCA_IDEAL_Counters(:,1,:);
    sinograms(:,sinogramColumn,:,3) = NoCSCA_DLR_Counters(:,1,:);
    sinograms(:,sinogramColumn,:,4) = NoCSCA_FDR_Counters(:,1,:);
    sinograms(:,sinogramColumn,:,5) = NoCSCA_PCS_Counters(:,1,:);
    sinograms(:,sinogramColumn,:,6) = NoCSCA_SR_Counters(:,1,:);
    
sinograms(:,sinogramColumn,:,7) = Add2x2Dy_Counters(:,1,:);
sinograms(:,sinogramColumn,:,8) = Add2x2St_Counters(:,1,:);
sinograms(:,sinogramColumn,:,9) = Add3x3Dy_Counters(:,1,:);
sinograms(:,sinogramColumn,:,10) = Add3x3St_Counters(:,1,:);
sinograms(:,sinogramColumn,:,11) = AddHybridDy_Counters(:,1,:);
sinograms(:,sinogramColumn,:,12) = AddHybridSt_Counters(:,1,:);

sinograms(:,sinogramColumn,:,13) = Sub2x2Dy_Counters(:,1,:);
sinograms(:,sinogramColumn,:,14) = Sub2x2St_Counters(:,1,:);
sinograms(:,sinogramColumn,:,15) = Sub3x3Dy_Counters(:,1,:);
sinograms(:,sinogramColumn,:,16) = Sub3x3St_Counters(:,1,:);
sinograms(:,sinogramColumn,:,17) = SubHybridDy_Counters(:,1,:);
sinograms(:,sinogramColumn,:,18) = SubHybridSt_Counters(:,1,:);

    end
end
%%
%clearvars -except sinograms
save(strcat(BASEFILENAME,'_',IDSTRING,'_','Sinogram_allOn_Thresholded_',num2str(NumOfThresh),'Thresh'),'sinograms','-v7.3')
strcat(BASEFILENAME,'_',IDSTRING,'_','Sinogram_allOn_Thresholded_',num2str(NumOfThresh),'Thresh')

%%
for ij = 1:NumOfApproaches
for jj = 1:NumOfThresh
    figure(ij)
mesh(sinograms(:,:,jj,ij))
view(0, 90)
title(strcat(num2str(ij), " ", num2str(jj)))
pause(1)
end
end
%%

for window = 1:NumOfThresh-1
    for Approach = 1:NumOfApproaches
        sinograms(:,:,window,Approach) = sinograms(:,:,window,Approach) - sinograms(:,:,window+1,Approach);
    end
end
sinograms(:,:,NumOfThresh,:) = sinograms(:,:,NumOfThresh,:);
%%
%clearvars -except sinograms
save(strcat(BASEFILENAME,'_',IDSTRING,'_','Sinogram_allOn_Windowed_',num2str(NumOfThresh),'Thresh'),'sinograms','-v7.3')
strcat(BASEFILENAME,'_',IDSTRING,'_','Sinogram_allOn_Windowed_',num2str(NumOfThresh),'Thresh')
%%
for ij = 1:NumOfApproaches
for jj = 1:NumOfThresh
    figure(ij)
mesh(sinograms(:,:,jj,ij))
view(0, 90)
title(strcat(num2str(ij), " ", num2str(jj)))
pause(1)
end
end
