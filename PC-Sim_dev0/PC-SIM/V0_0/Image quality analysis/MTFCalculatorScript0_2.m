%{
clc
clear
close all
% Ensure you have the 'npy-matlab' library to read numpy arrays in MATLAB
% You can get it from: https://github.com/kwikteam/npy-matlab
addpath('C:\Program Files\MATLAB\npy-matlab-master\npy-matlab-master\npy-matlab\');

% Load the image data
AllImages_data = readNPY('FBP_Images_array_FullPhantom_Smaller_RealisticContrast_10e7_allOn_1.npy');
image_data(:,:) = AllImages_data(1,1,:,:);
image_data = image_data - (min(min(image_data)));
image_data = image_data / (max(max(image_data)));
%}
%%
clear
clc
close all
% Prompt user to load the CT image
[file, path] = uigetfile({'*.png;*.jpg;*.tiff;*.dcm'}, 'Select the CT image');
ct_image_uint = imread(fullfile(path, file));

% Check if the image is color and convert to grayscale if necessary
if size(ct_image_uint, 3) == 3
    % Convert to grayscale
    ct_image_gray = rgb2gray(ct_image_uint);
else
    % Directly assign grayscale image
    ct_image_gray = ct_image_uint;
end

% Convert to double precision
ct_image = double(ct_image_gray);

image_data=ct_image;






%%
% Display the image and select the circular insert region
figure(1);
imshow(image_data, []);
imcontrast()
pause;
title('Select the circular insert region using an ellipse');
h = imellipse;
wait(h);
mask = createMask(h);

% Extract the region of the circular insert
phantom_region = image_data;
phantom_region(~mask) = 0;

% **Step 3: Normalize the image data**
% Prompting user to select ROIs for air and phantom
figure(2); imshow(image_data, []);
title('Select ROI for air');
air_roi = drawrectangle('Label','Air ROI');
wait(air_roi);
air_mask = createMask(air_roi);
x_air = mean(image_data(air_mask));

figure(3); imshow(image_data, []);
title('Select ROI for uniform phantom');
phantom_roi = drawrectangle('Label','Phantom ROI');
wait(phantom_roi);
phantom_mask = createMask(phantom_roi);
x_phantom = mean(image_data(phantom_mask));

% Apply the normalization formula
normalisedMTFCalculationSubRegion = (image_data - x_air) / (x_phantom - x_air);

% **Step 4: Identify the cylindrical phantom and find its center**
threshold_level = graythresh(phantom_region);
binary_image = imbinarize(phantom_region, threshold_level);

% Calculate the centroid of the circular region
props = regionprops(binary_image, 'Centroid', 'EquivDiameter');
center = props.Centroid;
radius = props.EquivDiameter / 2;



%Calculate ESF
    %Calculate centre point of circle

        viscircles(center, radius,'EdgeColor','r');
        max_radius = radius*1.75;
        viscircles(center, max_radius,'EdgeColor','b');
%%        

showPlots = 0;
    %Make radial lines to outside of ring
% Example: Assuming you have an image stored in 'imageMatrix'
angularStepSize = 0.1;
theta_degrees = 0:angularStepSize:360-angularStepSize; % Angles in degrees
theta_radians = deg2rad(theta_degrees); % Convert to radians

% Convert polar coordinates to Cartesian coordinates
SamplingFrequency = 1
intensityValues(1:360/angularStepSize,1:uint16(max_radius*SamplingFrequency)) = 0;
figure(4)
hold off
figure(5)
hold off
movegui("northwest")
figure(6)
hold off
movegui("southwest")



for rho = 0:1/SamplingFrequency:max_radius
    rho
[x, y] = pol2cart(theta_radians, rho);

%work out radius in pol coordinates
figure(4)
if showPlots == 1
imshow(normalisedMTFCalculationSubRegion)
hold on
viscircles(center, radius,'EdgeColor','r');
viscircles(center, max_radius,'EdgeColor','b');
viscircles(center, rho,'EdgeColor','g');
end
movegui("center")

% Sample intensity values
intensityValues(1:360/angularStepSize,uint16(rho*SamplingFrequency)+1) = interp2(normalisedMTFCalculationSubRegion, center(1) +x, center(2)+y, 'linear');
if showPlots == 1
figure(5)
plot(intensityValues(:,uint16(rho*SamplingFrequency)))
hold on
title(rho*SamplingFrequency)
figure(6)
polarplot(theta_radians, intensityValues);
title('Intensity Profile Along Radial Lines');
pause(0.01)
end

end

%%

figure(5)
movegui("northeast")
hold off
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

figure(6)
movegui("southeast")
hold off
figure(6)
plot(mean(intensityValues(:,:),1))

%%

%Differentiate ESF to get LSF
% this is the "finite difference" derivative. Note it is  one element shorter than y and x
clear x
x = [1:length(mean(intensityValues(:,:),1))];
yd = abs(diff(mean(intensityValues(:,:),1))./diff(x));
LSF = yd;
% this is to assign yd an abscissa midway between two subsequent x
 xd = (x(2:end)+x(1:(end-1)))/2;
 % this should be a rough plot of derivative
 figure(7)
 hold off
 plot(xd,yd)
hold on



%%
L = length(LSF); % Length of the LSF

%%If Hann filtering
hann_window = hann(L); % Generate the Hann window
LSF_hann = LSF .* hann_window'; % Apply the window
LSF = LSF_hann;
%%If Hann filtering

figure(7)
hold on
plot(xd,LSF_hann)
hold off

P2 = abs(fft(LSF / L)); % Compute the FFT
P1 = P2(1:floor(L/2)+1); % Take only the positive frequencies
%P1(2:end-1) = 2 * P1(2:end-1); % Double the amplitudes (except DC and Nyquist)

nP1 = P1 - min(P1); % Shift to non-negative values
nP1 = nP1 ./ max(nP1); % Normalize to [0, 1]

Fs = 2; % Sampling frequency (pixels/lp)
f = Fs * (0:(L/2)) / L; % Normalized spatial frequency

figure(8)
plot(f(1:end), nP1(1:end), 'o-','linewidth', 3);
title('Radial MTF');
xlabel('Normalized Spatial Frequency');
ylabel('MTF');


%% NEED TO FIGURE OUT WAY TO INCRREASE SAMPLING FREQUENCY IN THE INTERPOLATION OF THE POLAR COORDINATES SOME HOW PERHAPS?
%%




function [SubRegion] = areaSampler(image, label)
figure(1)
imshow(image)
title(label)
f=figure(1);
set(f,'Position',get(0,'screensize'))
sampledArea = getrect;
left_x = sampledArea(1);
bottom_y = sampledArea(2);
width_dx = sampledArea(3);
height_dy = sampledArea(4);
rectangle('Position', [left_x, bottom_y, width_dx, height_dy], 'EdgeColor', 'c', 'LineWidth', 2);
SubRegion = image(bottom_y:(bottom_y + height_dy), left_x:(left_x + width_dx));
close(1)
end

function bigPlot(image, label)
figure(2)
imshow(image)
title(label)
f=figure(2);
set(f,'Position',get(0,'screensize'))
end