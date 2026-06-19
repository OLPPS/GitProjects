clear; clc; close all;

%% =============================
% USER GUI INPUT
%% =============================

% --- CSA ---
prompt = {'Enter CSA capacitance Cf (pF):', ...
          'Enter CSA feedback resistor Rf (MOhm):'};
definput = {'1.4', '100'};  % realistic
answer = inputdlg(prompt,'CSA Parameters',[1 50], definput);

Cf = str2double(answer{1}) * 1e-12;
Rf = str2double(answer{2}) * 1e6;

% --- CR stage ---
prompt = {'Enter CR capacitance Cc (pF):', ...
          'Enter CR resistor Rc (kOhm):'};
definput = {'1', '10'};
answer = inputdlg(prompt,'CR Stage Parameters',[1 50], definput);

C_cr = str2double(answer{1}) * 1e-12;
R_cr = str2double(answer{2}) * 1e3;

% --- RC stages ---
prompt = {'Number of RC stages (n):'};
answer = inputdlg(prompt,'RC Stages',[1 50], {'2'});
n = str2double(answer{1});

C_rc = zeros(1,n);
R_rc = zeros(1,n);

if n > 0
    prompt = cell(1,2*n);
    definput = cell(1,2*n);

    for i=1:n
        prompt{2*i-1} = sprintf('RC%d Capacitance (pF):',i);
        prompt{2*i}   = sprintf('RC%d Resistance (kOhm):',i);
        definput{2*i-1} = '1';
        definput{2*i}   = '10';
    end

    answer = inputdlg(prompt,'RC Values',[1 50], definput);

    for i=1:n
        C_rc(i) = str2double(answer{2*i-1}) * 1e-12;
        R_rc(i) = str2double(answer{2*i})   * 1e3;
    end
end

%% =============================
% TIME SETUP
%% =============================

t_max = 200e-9;
dt    = 1e-9;
t     = 0:dt:t_max;

%% =============================
% LOAD MAP
%% =============================

data = load_preprocessed_series('HEXITEC_MHZ_CdTeFromAttemptedCZT_250micron_2mm_1000V_centre_both_%dns',35);
[Qmap_txyz, nx, ny, nz] = reshape_to_4D(data);

%% =============================
% TEST VOXEL
%% =============================

ix=round(nx/2); iy=round(ny/2); iz=round(nz/2);
Q = squeeze(Qmap_txyz(:,ix,iy,iz));

%% =============================
% CONVOLUTION MODEL
%% =============================

tic;
V_conv = convolution_chain(Q, t, Cf, Rf, C_cr, R_cr, C_rc, R_rc);
t_conv = toc;

%% =============================
% PCA ANALYTICAL MODEL
%% =============================

tic;
[V_analytic, pca_model] = pca_analytical_chain(Q, t, Cf, Rf, C_cr, R_cr, C_rc, R_rc);
t_analytic = toc;

%% =============================
% FIX BASELINE
%% =============================

V_conv = V_conv - V_conv(1);
%V_conv = -V_conv;

V_analytic = V_analytic - V_analytic(1);
%V_analytic = -V_analytic;

%% =============================
% COMPARISON PLOT
%% =============================

figure;
plot(t, V_conv,'b--','LineWidth',2); hold on;
plot(t, V_analytic,'r','LineWidth',2);
legend('Convolution','PCA Analytical');

xlabel('Time (s)');
ylabel('Signal');
title('Shaper Output Comparison');

disp(['Conv time: ',num2str(t_conv)]);
disp(['Analytic time: ',num2str(t_analytic)]);