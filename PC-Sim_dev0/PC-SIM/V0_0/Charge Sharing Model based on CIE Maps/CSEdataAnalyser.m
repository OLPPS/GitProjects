clear
close all
clc
data = dlmread('C:\Users\oliep\OneDrive\Work\PDTF 3 - 2023-2027 (EPSRC SPCS Project)\SPCS codes\edges_test_600um');
%data = dlmread('HEXITEC_MHZ_CdTeFromAttemptedCZT_250micron_2mm_1000V_centre_both_35ns');
%
% Extract the X, Y, Z coordinates
X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

% Find the range of X, Y, and Z
xRange = [min(X), max(X)];
yRange = [min(Y), max(Y)];
zRange = [min(Z), max(Z)];

fprintf('Range of X: [%f, %f]\n', xRange(1), xRange(2));
fprintf('Range of Y: [%f, %f]\n', yRange(1), yRange(2));
fprintf('Range of Z: [%f, %f]\n', zRange(1), zRange(2));

% Ask the user for the desired X value and time t
desiredX = input('Enter the desired constant X value: ');
desiredT = input('Enter the time step t (0 to 35): ');

% Validate the time step
if desiredT < 0 || desiredT > 35
    error('Invalid time step. Please enter a value between 0 and 35.');
end

% Calculate the column index for the specified time t (4th column onward)
valueColumn = 4 + desiredT;

% Extract the slice where X equals the desired value
slice = data(abs(X - desiredX) < 1e-6, :); % Tolerance for floating-point comparison

% Extract the Y, Z, and value data for the given time step
Y_slice = slice(:, 2);
Z_slice = slice(:, 3);
values = slice(:, valueColumn);

% Display the slice data
fprintf('Extracted slice at X = %f and t = %dns:\n', desiredX, desiredT);
result = [Y_slice, Z_slice, values];
%disp(result);

%
%result = [0,0,0; 0,1,1; 0,2,0; 1,0,1; 1,1,2; 1,2,1; 2,0,0; 2,1,1; 2,2,0;];
figure(1)
surf(result)
view(0,90)

clc
Y = result(:, 1); % First column for Y
Z = result(:, 2); % Second column for Z
Value = result(:, 3); % Third column for the values

Y_unique = unique(Y);
Z_unique = unique(Z);
[Y_grid, Z_grid] = meshgrid(Z_unique, Y_unique);

Value_grid = reshape(Value, length(Y_unique), length(Z_unique));

figure(2)
surf(Y_grid,Z_grid,Value_grid)
view(90,0)
pause(0.5)
%%
figure(3)
maxSignal = 0.93935;%0.95981;
for ii = 1:7:889
    acrossY(ii,:) = result(ii:ii+6,3);
plot(Y_unique,acrossY(ii,:))
hold on
for ij = 1:7
    totalmirroredY(ii,ij) = acrossY(ii,ij) + acrossY(ii,8-ij);
CSmirroredY(ii,ij) = acrossY(ii,8-ij);
    lostMethodY(ii,ij) = maxSignal - acrossY(ii,ij);
    totalLostMethodY(ii,ij) = maxSignal - acrossY(ii,ij) + acrossY(ii,ij);
end
plot(Y_unique,CSmirroredY(ii,:))
plot(Y_unique,lostMethodY(ii,:))
plot(Y_unique,totalmirroredY(ii,:))
plot(Y_unique,totalLostMethodY(ii,:))
title(ii)
legend
hold off
pause(0.1)
end


figure(4)
maxSignal = 0.93935;%0.95981;
for ii = 1:7
    acrossZ(ii,:) = result(ii:7:889,3);
plot(Z_unique,acrossZ(ii,:))
hold on
for ij = 1:127
    totalmirroredZ(ii,ij) = acrossZ(ii,ij) + acrossZ(ii,128-ij);
CSmirroredZ(ii,ij) = acrossZ(ii,128-ij);
    lostMethodZ(ii,ij) = maxSignal - acrossZ(ii,ij);
    totalLostMethodZ(ii,ij) = maxSignal - acrossZ(ii,ij) + acrossZ(ii,ij);
end
plot(Z_unique,CSmirroredZ(ii,:))
plot(Z_unique,lostMethodZ(ii,:))
plot(Z_unique,totalmirroredZ(ii,:))
plot(Z_unique,totalLostMethodZ(ii,:))
title(ii)
legend
hold off
pause(1)
end

%%
result(:,1)

%%
% Example input array
step = 127; % Step size
rows = size(result, 1); % Total number of rows

% Preallocate the reshaped array
reshapedArray = [];

% Loop through the starting rows (1 to step)
for startRow = 1:step
    % Extract rows starting at 'startRow' with a step of 'step'
    extractedRows = result(startRow:step:rows, :);
    % Append the extracted rows to the reshaped array
    reshapedArray = [reshapedArray; extractedRows];
end




%%
A = unique(data(:,1));
size(A)

%%
clc
Y = testresult(:, 1); % First column for Y
Z = testresult(:, 2); % Second column for Z
Value = testresult(:, 3); % Third column for the values

Y_unique = unique(Y);
Z_unique = unique(Z);
[Y_grid, Z_grid] = meshgrid(Y_unique, Z_unique);

Value_grid = reshape(Value, length(Z_unique), length(Y_unique));

C = zeros(size(Value_grid) - [1, 1]); % Initialize matrix for averages
for i = 1:size(C, 1)
    for j = 1:size(C, 2)
        % Average of the four surrounding points
        C(i, j) = mean([Value_grid(i, j), Value_grid(i+1, j), Value_grid(i, j+1), Value_grid(i+1, j+1)]);
    end
end

% Adjust X and Y to match the size of C
Y = X(1:end-1, 1:end-1);
Z = Y(1:end-1, 1:end-1);

% Plot the surface with the new color matrix
figure(1);
surf(Y_grid(1:end-1,end-1), Z_grid(1:end-1,end-1), Value_grid(1:end-1, 1:end-1), C);
%surf(Y_grid,Z_grid,Value_grid)
colorbar;
xlabel('Y_grid');
ylabel('Z_grid');
zlabel('Value_grid');
title('Surface Plot with Averaged Colors');
view(0,90)
%}



