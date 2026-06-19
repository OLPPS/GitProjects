clear; clc; close all;

%% =============================
% USER GUI INPUT
%% =============================

prompt = {'This model assumes a CR-RC^n shaper. Please enter an integer value for n:'};
dlgtitle = 'Shaper Configuration';
dims = [1 50];
definput = {'2'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
n = str2double(answer{1});

% Now request C and R values
prompt = cell(1,2*n);
for i = 1:n
    prompt{2*i-1} = ['Enter C for stage ' num2str(i) ' (pF):'];
    prompt{2*i}   = ['Enter R for stage ' num2str(i) ' (kOhm):'];
end

dlgtitle = 'Enter RC parameters [DEFAULTS: C = 1, R = 1]';
dims = [1 50];
definput = repmat({'1'},1,2*n);
answer = inputdlg(prompt, dlgtitle, dims, definput);

Cvals = zeros(1,n);
Rvals = zeros(1,n);

for i = 1:n
    Cvals(i) = str2double(answer{2*i-1})*10^(-12);
    Rvals(i) = str2double(answer{2*i})*10^3;
end

%% =============================
% FIXED PARAMETERS (EDIT HERE)
%% =============================

t_max = 200e-9;     % waveform duration (s)
dt    = 1e-9;       % timestep
t     = 0:dt:t_max;

nt = length(t);

% Example 4D Prettyman map (REPLACE WITH YOUR DATA)
%nx=30; ny=30; nz=50;
%Qmap = rand(nx,ny,nz,nt);  % placeholder

%data = dlmread('C:\Users\oliep\OneDrive\Work\PDTF 3 - 2023-2027 (EPSRC SPCS Project)\SPCS codes\centre_test_600um');

data = load_preprocessed_series('HEXITEC_MHZ_CdTeFromAttemptedCZT_250micron_2mm_1000V_centre_both_%dns',35);


[Qmap_txyz, nx, ny, nz] = reshape_to_4D(data);

nt_CIEmap = size(Qmap_txyz,1)


%% =============================
% RESHAPE MAP
%% =============================

%Qmap_txyz = reshape_to_txyz(Qmap);

%% =============================
% VISUALISATION LOOP
%% =============================

figure(5);

mid_x = round(nx/2);
mid_y = round(ny/2);
mid_z = round(nz/2);
%{
for ti = 1:nt_CIEmap

    clf;

    % XY slice
    subplot(2,2,1)
    imagesc(squeeze(Qmap_txyz(ti,:, :, mid_z)));
    title('XY slice')
    colorbar

    % YZ slice
    subplot(2,2,2)
    imagesc(squeeze(Qmap_txyz(ti,mid_x, :, :)));
    title('YZ slice')
    colorbar

    % XZ slice
    subplot(2,2,3)
    imagesc(squeeze(Qmap_txyz(ti,:, mid_y, :)));
    title('XZ slice')
    colorbar

    % Time info
    subplot(2,2,4)
    axis off
    text(0.3,0.5, sprintf('Time step: %d\nTime = %.2e s',ti,t(ti)), ...
        'FontSize',12)

    drawnow;
    pause(0.2);

end
%%
figure(1)
for zz = 1:nz

    clf;
% XY slice
    imagesc(squeeze(Qmap_txyz(end,:, :, zz)));
    title('XY slice')
    colorbar
    pause(0.2)
end

figure(2)
for xx = 1:nx

    clf;
% YZ slice
    imagesc(squeeze(Qmap_txyz(end,xx, :, :)));
    title('XY slice')
    colorbar
    pause(0.2)
end

figure(3)
for yy = 1:ny

    clf;
% XZ slice
    imagesc(squeeze(Qmap_txyz(end,:, yy,:)));
    title('XY slice')
    colorbar
    pause(0.2)
end
%}
%% =============================
% TEST SINGLE VOXEL
%% =============================

ix=mid_x; iy=mid_y; iz=mid_z;

Q = squeeze(Qmap_txyz(:,ix,iy,iz));

%% ---- Analytical approach ----
tic;
V_analytic = analytical_chain(Q, t, Rvals, Cvals);
t_analytic = toc;

%% ---- Convolution approach ----
tic;
V_conv = convolution_chain(Q, t, Rvals, Cvals);
t_conv = toc;

%% =============================
% PLOT COMPARISON
%% =============================

V_conv = V_conv - V_conv(1);
V_conv = - V_conv;

figure;
plot(t, V_analytic, 'r', 'LineWidth',2); hold on;
plot(t, V_conv, '--b', 'LineWidth',2);
legend('Analytical','Convolution');
xlabel('Time (s)');
ylabel('Signal');
title('Waveform Comparison');

disp(['Analytical time: ' num2str(t_analytic)]);
disp(['Convolution time: ' num2str(t_conv)]);



function [V4, Nx, Ny, Nz] = reshape_to_4D(data)
% data is MxN: [x, y, z, v1, v2, ..., vT]

    % Extract coordinates
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);

    [x1,x2] = maxmin(x)
    [y1,y2] = maxmin(y)
    [z1,z2] = maxmin(z)

    % Extract values
    V = data(:,4:end);   % M x T
    T = size(V,2);

    % Unique coordinate grids
    ux = unique(x);
    uy = unique(y);
    uz = unique(z);

    Nx = numel(ux);
    Ny = numel(uy);
    Nz = numel(uz);

    % Preallocate output
    V4 = nan(T, Nx, Ny, Nz);

    % Build index maps
    [~, ix] = ismember(x, ux);
    [~, iy] = ismember(y, uy);
    [~, iz] = ismember(z, uz);

    % Fill the 4D array
    for row = 1:length(x)
        V4(:, ix(row), iy(row), iz(row)) = V(row, :).';
    end
end
