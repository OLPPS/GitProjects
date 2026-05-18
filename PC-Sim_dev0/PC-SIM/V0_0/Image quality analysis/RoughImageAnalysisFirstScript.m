% Prompt user to load the CT image
[file, path] = uigetfile({'*.png;*.jpg;*.tiff;*.dcm'}, 'Select the CT image');
ct_image_uint = imread(fullfile(path, file));

% Check if the image is color and convert to grayscale if necessary
if size(ct_image_uint, 3) == 3
    % Define ct_image_color and convert to grayscale
%    ct_image_color = ct_image_uint; %REMOVE
    ct_image_gray = rgb2gray(ct_image_uint);
else
    % Directly assign grayscale image
    ct_image_gray = ct_image_uint;
end

% Convert to double precision
ct_image = double(ct_image_gray);
%%
% Display the image with initial windowing
h = imshow(ct_image_gray, []);
title('Adjust Contrast Using the Tool');

% Open the imcontrast tool for manual adjustment
imcontrast(h);

% Wait for the user to adjust the contrast
disp('Adjust the contrast using the tool, then press Enter to continue...');
pause;

% Retrieve the adjusted contrast limits from the imcontrast tool
%cmap = get(h, 'CDataMapping'); %REMOVE
clim = get(gca, 'CLim'); % Get the current contrast limits (CLim) from the axes

% Apply the contrast limits to the original image data
scaled_ct_image = (ct_image - clim(1)) / (clim(2) - clim(1));
scaled_ct_image(scaled_ct_image < 0) = 0;
scaled_ct_image(scaled_ct_image > 1) = 1;

% Display the contrast-adjusted image
figure, imshow(scaled_ct_image), title('Contrast-Adjusted CT Image');


%%
% Step 1: Calculate Contrast-to-Noise Ratio (CNR)
disp('Step 1: Please select two regions - one with the structure of interest and one background region.');
roi_structure = drawrectangle('Label', 'Structure of Interest', 'Color', 'r');
roi_background = drawrectangle('Label', 'Background Region', 'Color', 'b');

% Map the selected regions to the original image
structure_values = ct_image(round(roi_structure.Position(2)):(round(roi_structure.Position(2))+round(roi_structure.Position(4))), ...
                            round(roi_structure.Position(1)):(round(roi_structure.Position(1))+round(roi_structure.Position(3))));
background_values = ct_image(round(roi_background.Position(2)):(round(roi_background.Position(2))+round(roi_background.Position(4))), ...
                             round(roi_background.Position(1)):(round(roi_background.Position(1))+round(roi_background.Position(3))));

% Calculate mean and standard deviation for the regions
mean_structure = mean(structure_values(:));
mean_background = mean(background_values(:));
std_background = std(background_values(:));

% Check if the selected background region has a standard deviation of zero
if std_background == 0
    error('The selected background region has a standard deviation of zero. Please select a different background region.');
end

% Calculate CNR
CNR = (mean_structure - mean_background) / std_background;
fprintf('Contrast-to-Noise Ratio (CNR): %.2f\n', CNR);

% Step 2: Calculate Signal-to-Noise Ratio (SNR)
disp('Step 2: Please select a region with the signal of interest.');
roi_signal = drawrectangle('Label', 'Signal Region', 'Color', 'g');

% Map the selected region to the original image
signal_values = ct_image(round(roi_signal.Position(2)):(round(roi_signal.Position(2))+round(roi_signal.Position(4))), ...
                         round(roi_signal.Position(1)):(round(roi_signal.Position(1))+round(roi_signal.Position(3))));

% Calculate mean and standard deviation for the signal region
mean_signal = mean(signal_values(:));
std_signal = std(signal_values(:));

% Check if the selected signal region has a standard deviation of zero
if std_signal == 0
    error('The selected signal region has a standard deviation of zero. Please select a different signal region.');
end

% Calculate SNR
SNR = mean_signal / std_signal;
fprintf('Signal-to-Noise Ratio (SNR): %.2f\n', SNR);

% Step 3: Calculate Entropy
Entropy = entropy(ct_image);
fprintf('Entropy: %.2f\n', Entropy);

% Step 4: Calculate Local Variance
disp('Step 4: Calculating Local Variance.');
local_var = stdfilt(ct_image).^2;
mean_local_var = mean(local_var(:));
fprintf('Local Variance: %.2f\n', mean_local_var);

% Step 5: Calculate Noise Power Spectrum (NPS)
disp('Step 5: Calculating Noise Power Spectrum (NPS).');
fft_image = fftshift(fft2(ct_image));
nps = abs(fft_image).^2;
mean_nps = mean(nps(:));
fprintf('Noise Power Spectrum (NPS): %.2f\n', mean_nps);

% Step 6: Calculate Gradient Magnitude
disp('Step 6: Calculating Gradient Magnitude.');
[Gx, Gy] = gradient(ct_image);
gradient_magnitude = sqrt(Gx.^2 + Gy.^2);
mean_gradient_magnitude = mean(gradient_magnitude(:));
fprintf('Gradient Magnitude: %.2f\n', mean_gradient_magnitude);

% Step 7: Calculate Uniformity Index
disp('Step 7: Calculating Uniformity Index.');
uniformity_index = std(ct_image(:)) / mean(ct_image(:));
fprintf('Uniformity Index: %.2f\n', uniformity_index);

% Step 8: Calculate SSIM based on Gaussian-blurred image
disp('Step 8: Calculating SSIM based on Gaussian-blurred image.');
gaussian_blurred_image = imgaussfilt(ct_image, 2); % Apply Gaussian blur with a sigma of 2
SSIM_val = ssim(ct_image, gaussian_blurred_image);
fprintf('Structural Similarity Index (SSIM): %.2f\n', SSIM_val);

% Output the results to the user
disp('--------------------');
disp('Results:');
fprintf('Contrast-to-Noise Ratio (CNR): %.2f\n', CNR);
fprintf('Signal-to-Noise Ratio (SNR): %.2f\n', SNR);
fprintf('Entropy: %.2f\n', Entropy);
fprintf('Local Variance: %.2f\n', mean_local_var);
fprintf('Noise Power Spectrum (NPS): %.2f\n', mean_nps);
fprintf('Gradient Magnitude: %.2f\n', mean_gradient_magnitude);
fprintf('Uniformity Index: %.2f\n', uniformity_index);
fprintf('Structural Similarity Index (SSIM): %.2f\n', SSIM_val);
disp('--------------------');

% Inform the user that the analysis is complete
disp('Analysis complete. The results are displayed above.');
