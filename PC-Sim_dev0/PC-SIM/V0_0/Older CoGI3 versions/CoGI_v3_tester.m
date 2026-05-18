clear
clc

MapAddresses{1} = 'thickness_1pt5mm_200MicronPitch_Central_Max_saved';
MapAddresses{2} = 'thickness_1pt5mm_200MicronPitch_Central_Max_saved';
MapAddresses{3} = 'thickness_1pt5mm_200MicronPitch_Central_Max_saved';
MapAddresses{4} = 'thickness_1pt5mm_200MicronPitch_LeftEdge_Max_saved';
MapAddresses{5} = 'thickness_1pt5mm_200MicronPitch_TopLeftCorner_Max_saved';


% Define kernel size and standard deviation
kernelSize = 11;       % Must be odd for symmetry
sigma = (kernelSize - 1)/6;             % Standard deviation

% Create coordinate grid centered at zero
r = (kernelSize - 1) / 2;
[x, y, z] = ndgrid(-r:r, -r:r, -r:r);

% Spherical Gaussian formula
G = exp(-(x.^2 + y.^2 + z.^2) / (2 * sigma^2));

% Normalize kernel to preserve energy
G = G / sum(G(:));

orientation = 'stdrev';

if orientation == 'stdrev'
for CSMap = 1:9
    %Load in map data and rotate as needed
    switch CSMap
        case 1
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 0+1), [3 1 2]);
        case 2
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 0), [3 1 2]);
        case 3
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 1+3), [3 1 2]);
        case 4
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 3+2), [3 1 2]);
        case 5
            CIEMapAddress = MapAddresses{1};
            load(CIEMapAddress);
        case 6
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 1+2), [3 1 2]);
        case 7
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 3+3), [3 1 2]);
        case 8
            CIEMapAddress = MapAddresses{4};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 2), [3 1 2]);
        case 9
            CIEMapAddress = MapAddresses{5};
            load(CIEMapAddress);
            CIEMap = permute(rot90(permute(CIEMap, [2 3 1]), 2+1), [3 1 2]);
    end
   storedMaps{CSMap} = CIEMap;
end
end




figure(333)

for depth = 103:103
    tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'none');
    for k = 1:9
    % Create subplot in a 3x3 grid
nexttile
    %subplot(3,3,k);

    A = storedMaps{k};
X = squeeze(A(depth,:,:));
%X = squeeze(sum(A,1));
surf(X)
imagesc(X,[-0.1,1.1]);
caxis([-0.1 1.1]);
shading interp;
view(0,90)
%    title(['Plot ' num2str(k) ' ' num2str(depth)]);
axis off
end
pause(0.1)
colormap('turbo');
end

%%
for depth = 1:201
A = storedMaps{5};
figure(111)
X = squeeze(A(depth,:,:));
%X = squeeze(sum(A,1));
surf(X);
shading interp;
view(0,90)
title(depth)
pause(0.1)
end


%%
% ParamsFile holds various parameters needed at this stage, e.g. which CSCAs to use.
% GeometryFile defines the geometry of the detector with respect to the x-ray source in STANDARDISED format.
% NoiseParametersFile imports information regarding the noise profile of the detector.
% MCDataMatFileName contains the MC data in STANDARDISED format.
% t1, and t2 are the start and end times of the slice being processed.
% RANDOMSEED is set to allow reproducibillity of the simulation and fair comparisson between various CSCAs or ACTS etc.
% CIEtolerance is a threshold which can be set such that events with a CIE value below this point are not processed further.
% SectionsToRun this is a vector which stores which sections of CoGI to run. See documentation for more details.
% MapAddresses is an array which stores the addresses of the CIE maps for centre, edge, corner, CSEdge and CSCorner in STANDARISED format, in that order.
% VERBOSE sets text output (1 = a lot, default = basic)
% DEBUGGING sets whether debugging mode is on (see documentation for more details.)
% Label is the base file name for outputs.





clear
clc
MapAddresses{1} = 'thickness_1pt5mm_200MicronPitch_Central_Max_saved';
MapAddresses{2} = 'thickness_1pt5mm_200MicronPitch_Central_Max_saved';
MapAddresses{3} = 'thickness_1pt5mm_200MicronPitch_Central_Max_saved';
MapAddresses{4} = 'thickness_1pt5mm_200MicronPitch_LeftEdge_Max_saved';
MapAddresses{5} = 'thickness_1pt5mm_200MicronPitch_TopLeftCorner_Max_saved';

tic
Start = datetime;

for ij =  1:1


%CoGI_v3_1_00_CERN_RunOne('PARAMS.mat',"GEOMETRIESFile_W0pt17_P0pt2_D1pt5_Xo0_Yo0_Zo50_MX36pt25_MY0pt3_Xpix362_Ypix3_rev","NoiseParameters_CERN_RunOne",'WO_Sm_RC_Meta_10e6_0_10e6__FormattedCorrectlyAndCompressed.mat',0,1/2000,12345,0.001, [1 1 1], MapAddresses,-1, 1,1,strcat('TEST1_',num2str(ij)),[2:2:10]);
CoGI_v3_1_00('PARAMS.mat',"GEOMETRIESFile_W0pt17_P0pt2_D1pt5_Xo0_Yo0_Zo50_MX36pt25_MY0pt3_Xpix362_Ypix3_rev","NoiseParameters_noisey",'WO_Sm_RC_Meta_10e6_0_10e6__FormattedCorrectlyAndCompressed.mat',0,1/2000,12345,0.001, [1 1 1], MapAddresses,-1, 1,1,strcat('NoiseTest1_',num2str(ij)),[2:1:20]);

end
Stop = datetime;
Stop - Start
toc
%%
 % Create matfile object
    mfile = matfile('WO_Sm_RC_Meta_10e6_0_10e6__FormattedCorrectlyAndCompressed.mat');
%%
 % Access target column
    times = mfile.('EventsData')(:, 4);  % Loads only column number TimeColNumber into RAM

    %%
    clc
    PauseTime = 1
    figure(1)
    hold off
    for ij = 1:length(Thresholds)
    plot(ArrayOfCountsToAdd(:,2,ij),'DisplayName',strcat(num2str(ij*2),' keV threshold'))
    xlim([0,700])
    ylim([0 15000])
    xlabel('Pixel ID Number')
    ylabel('Dark field counts')
    hold on
    hold off
legend
legendHandle = legend; % Get the legend handle
legendHandle.Location = 'north'
legendHandle.FontSize = 16; % Set the font size to 14
    pause(PauseTime)
    end
%%
figure(2)
    hold off
    for ij = 1:1
    surf(ArrayOfCountsToAdd(95:133,:,ij),'DisplayName','STD')
    hold on
legend
view(0,90)
    pause(PauseTime)
    end

    %%
    for i = 1:1
    plot(NoCSCA_SR_Counters(:,2,ij),'DisplayName','SR')
legend
    pause(PauseTime)
    plot(NoCSCA_PCS_Counters(:,2,ij),'DisplayName','PCS')
legend
    pause(PauseTime)
    plot(NoCSCA_IDEAL_Counters(:,2,ij),'DisplayName','IDEAL')
legend
    pause(PauseTime)
    plot(NoCSCA_FDR_Counters(:,2,ij),'DisplayName','FDR')
legend
    pause(PauseTime)
    plot(NoCSCA_DLR_Counters(:,2,ij),'DisplayName','DLR')
legend
    pause(PauseTime)
    plot(SubHybridSt_Counters(:,2,ij),'DisplayName','-HySt')
legend
    pause(PauseTime)
    plot(SubHybridDy_Counters(:,2,ij),'DisplayName','-HyDy')
legend
    pause(PauseTime)
    plot(AddHybridSt_Counters(:,2,ij),'DisplayName','+HySt')
legend
    pause(PauseTime)
    plot(AddHybridDy_Counters(:,2,ij),'DisplayName','+HyDy')
legend
    pause(PauseTime)
    plot(Sub2x2St_Counters(:,2,ij),'DisplayName','-2St')
legend
    pause(PauseTime)
    plot(Sub2x2Dy_Counters(:,2,ij),'DisplayName','-2Dy')
legend
    pause(PauseTime)
    plot(Add2x2St_Counters(:,2,ij),'DisplayName','+2St')
legend
    pause(PauseTime)
    plot(Add2x2Dy_Counters(:,2,ij),'DisplayName','+2Dy')
legend
    pause(PauseTime)
    plot(Sub3x3St_Counters(:,2,ij),'DisplayName','-3St')
legend
    pause(PauseTime)
    plot(Sub3x3Dy_Counters(:,2,ij),'DisplayName','-3Dy')
legend
    pause(PauseTime)
    plot(Add3x3St_Counters(:,2,ij),'DisplayName','+3St')
legend
    pause(PauseTime)
    plot(Add3x3Dy_Counters(:,2,ij),'DisplayName','+3Dy')
hold off    
    end

    %%
    figure(2)
    hold off
    plot([10:120],squeeze(sum(sum(NoCSCA_STD_Counters(:,:,:),2),1)),'DisplayName','STD')
    hold on
    plot([10:120],squeeze(sum(sum(NoCSCA_SR_Counters(:,:,:),2),1)),'DisplayName','SR')
    plot([10:120],squeeze(sum(sum(AddHybridDy_Counters(:,:,:),2),1)),'DisplayName','+HyDy')
        plot([10:120],squeeze(sum(sum(Sub3x3St_Counters(:,:,:),2),1)),'DisplayName','-3St')
            plot([10:120],squeeze(sum(sum(Add2x2St_Counters(:,:,:),2),1)),'DisplayName','+2St')
            legend