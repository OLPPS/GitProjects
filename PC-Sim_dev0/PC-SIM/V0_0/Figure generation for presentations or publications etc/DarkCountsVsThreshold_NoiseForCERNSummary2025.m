clear
load('NoiseTest1_1_WholeSimFinished_1S1_2S1_3S1v7pt3.mat')

%{
figure(17)
hold off
plot((NoCSCA_STD_Counters(:,1,1)))
hold on
for ii = 2:5
pause(1)
plot((NoCSCA_STD_Counters(:,1,ii)))
end
%}

C(1) = 'r';
C(3) = 'b';
C(5) = 'k';
C(7) = 'g';
C(9) = 'm';

L{1} = '2 keV';
L{3} = '4 keV';
L{5} = '6 keV';
L{7} = '8 keV';
L{9} = '10 keV';

figure(18)
tiledlayout(1,5)
ax1 = nexttile([1 3])
hold off
%ArrayOfCountsToAdd(ArrayOfCountsToAdd(:,1,:) < 1) = 1;
data = ArrayOfCountsToAdd(:,:,1);
plot(ax1,data(:),'-','Color',C(1), 'LineWidth',2,'DisplayName',L{1})
hold on
for ii = 3:2:9
%pause(1)
data = ArrayOfCountsToAdd(:,:,ii);
plot(ax1,data(:),'-','Color',C(ii), 'LineWidth',2,'DisplayName',L{ii})
end
legend('FontSize',15,'NumColumns',5)
for jj = 1:19
   S(jj) = mean(ArrayOfCountsToAdd(:,1,jj));
end
xlabel('Pixel Number','FontSize',23)
ylabel('Dark Counts','FontSize',23)
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 1050])

ax2 = nexttile([1 2])
%nexttile(3)
semilogy(ax2,[2:2:10],S(1:2:9),'--','LineWidth',2)
xlabel('Threshold (keV)','FontSize',23)
ylabel('Dark Counts','FontSize',23)
set(gca, 'FontSize', 14, 'FontWeight', 'bold');


%%
% MATLAB script to plot Gaussian area above threshold N
% As N is swept from mean towards one end

%clear; clc;
%close all;

% Parameters of Gaussian
mu = 0;          % mean
sigma = 1;       % standard deviation

% Define Gaussian function
gauss = @(x) (1/(sigma*sqrt(2*pi))) * exp(-(x-mu).^2/(2*sigma^2));

% Range for plotting
x = linspace(mu-3*sigma, mu+3*sigma, 1000);

% Sweep threshold N from mean to +infinity
N_values = linspace(mu, mu+3*sigma, 1000);

% Preallocate area array
areas = zeros(size(N_values));

% Compute area above each threshold
for k = 1:length(N_values)
    N = N_values(k);
    % Integral of Gaussian from N to infinity
    areas(k) = integral(gauss, N, inf);
end

% Plot Gaussian curve with one example threshold
figure(1);
plot(x, gauss(x), '-ob', 'LineWidth', 1.5); hold on;
N_example = mu + sigma;
area_x = linspace(N_example, mu+3*sigma, 9);
area_y = gauss(area_x);
fill([area_x fliplr(area_x)], [area_y zeros(size(area_y))], ...
     'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('x'); ylabel('Gaussian PDF');
title('Gaussian with shaded area above threshold N');
legend('Gaussian PDF','Area above N');

% Plot area vs threshold
figure(2);
semilogy(9*N_values(112:end)/max(N_values(112:end)),areas(112:end)/max(areas(112:end)), '-r', 'LineWidth', 2, 'DisplayName','Simple gaussian model');
xlabel('Threshold N');
ylabel('Area above N');
grid on;
hold on
semilogy(S(1:9)/max(S(1:9)),'-ok','LineWidth',2, 'DisplayName','Modelled Noise')
xlabel('Threshold (keV)','FontSize',23)
ylabel('Normalised Dark Counts','FontSize',23)
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend

%%
% MATLAB script to plot Gaussian area above threshold N
% As N is swept from mean towards one end

%clear; clc;
close all;

% Parameters of Gaussian
mu = 0;          % mean
sigma = 1;       % standard deviation

% Define Gaussian function
gauss = @(x) (1/(sigma*sqrt(2*pi))) * exp(-(x-mu).^2/(2*sigma^2));

% Range for plotting
x = linspace(mu-3*sigma, mu+3*sigma, 1000);

% Sweep threshold N from mean to +infinity
N_values = linspace(mu, mu+3*sigma, 1000);

% Preallocate area array
areas = zeros(size(N_values));

% Compute area above each threshold
for k = 1:length(N_values)
    N = N_values(k);
    % Integral of Gaussian from N to infinity
    areas(k) = integral(gauss, N, inf);
end

figure(19)
tiledlayout(1,5)
ax1 = nexttile([1 3])
hold off
%ArrayOfCountsToAdd(ArrayOfCountsToAdd(:,1,:) < 1) = 1;
data = ArrayOfCountsToAdd(:,:,1);
plot(ax1,data(:),'-','Color',C(1), 'LineWidth',2,'DisplayName',L{1})
hold on
for ii = 3:2:9
%pause(1)
data = ArrayOfCountsToAdd(:,:,ii);
plot(ax1,data(:),'-','Color',C(ii), 'LineWidth',2,'DisplayName',L{ii})
end
legend('FontSize',15,'NumColumns',5)
for jj = 1:19
   S(jj) = mean(ArrayOfCountsToAdd(:,1,jj));
end
xlabel('Pixel Number','FontSize',23)
ylabel('Dark Counts','FontSize',23)
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 1050])

% Plot area vs threshold
ax2 = nexttile([1 2]);
semilogy(9*N_values(112:end)/max(N_values(112:end)),areas(112:end)/max(areas(112:end)), '-r', 'LineWidth', 2, 'DisplayName','Simple gaussian model');
xlabel('Threshold N');
ylabel('Area above N');
grid on;
hold on
semilogy(S(1:9)/max(S(1:9)),'-ok','LineWidth',2, 'DisplayName','Modelled Noise')
xlabel('Threshold (keV)','FontSize',23)
ylabel('Normalised Dark Counts','FontSize',23)
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend