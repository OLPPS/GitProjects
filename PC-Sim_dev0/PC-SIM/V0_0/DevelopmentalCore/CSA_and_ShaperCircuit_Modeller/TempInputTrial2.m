clearvars -except tcsa_exp tcsa_tri texp11 texp12 ttri11 ttri12 rcsa_exp rcsa_tri rexp11 rexp12 rtri11 rtri12; clc; close all;

%% =============================
% STAGE 1: STRUCTURE
%% =============================

% CSA type
choice = questdlg('Select CSA Type:', ...
    'CSA', ...
    'Exponential','Triangular','Exponential');

csa_mode = strcmp(choice,'Exponential');

% shaping
choice = questdlg('Apply shaping?', ...
    'Shaping','Yes','No','Yes');
use_shaper = strcmp(choice,'Yes');

n_CR = 0; n_RC = 0;

if use_shaper
    answ = inputdlg({'Number of CR stages','Number of RC stages'}, ...
        'Shaper structure',[1 50],{'1','1'});
    n_CR = str2double(answ{1});
    
    n_RC = str2double(answ{2});
end



% =============================
% STAGE 2: PARAMETER TYPES (UPDATED)
% =============================

labels = {};
types  = {};   % store 'tau' or 'rc' or 'tri'

% --- CSA ---
labels{end+1} = 'CSA';

% --- CR ---
for k = 1:n_CR
    labels{end+1} = sprintf('CR%d',k);
end

% --- RC ---
for k = 1:n_RC
    labels{end+1} = sprintf('RC%d',k);
end

%% =============================
% GLOBAL PARAMETER MODE
%% =============================

choice = questdlg( ...
    'How would you like to define CSA and shaper discharge behaviour?', ...
    'Parameter Mode', ...
    'All tau','All R & C','Mixed','All tau');

%% =============================
% ASSIGN TYPES
%% =============================

if strcmp(choice,'All tau')

    for i = 1:length(labels)
        if i==1 && ~csa_mode
            types{i} = 'tri';
        else
            types{i} = 'tau';
        end
    end

elseif strcmp(choice,'All R & C')

    for i = 1:length(labels)
        if i==1 && ~csa_mode
            types{i} = 'tri';
        else
            types{i} = 'rc';
        end
    end

else
    % ===== MIXED MODE =====

    prompt = cell(size(labels));

    for i=1:length(labels)
        if i==1 && ~csa_mode
            prompt{i} = sprintf('%s: type (tri)',labels{i});
        else
            prompt{i} = sprintf('%s: type (tau / rc)',labels{i});
        end
    end

    def = repmat({'tau'}, size(prompt));

    answ = inputdlg(prompt,'Parameter Types',[1 60],def);

    for i=1:length(labels)
        if i==1 && ~csa_mode
            types{i} = 'tri';
        else
            types{i} = lower(strtrim(answ{i}));
        end
    end

end


% =============================
% STAGE 3: PARAMETER VALUES
% =============================

prompt = {};
def = {};

for i=1:length(labels)

    name = labels{i};
    type = types{i};

    if strcmp(type,'tau')

        prompt{end+1} = sprintf('%s tau (ns)',name);
        def{end+1} = '20';

    elseif strcmp(type,'rc')

        prompt{end+1} = sprintf('%s R (kOhm)',name);
        def{end+1} = '100';

        prompt{end+1} = sprintf('%s C (pF)',name);
        def{end+1} = '0.2';

    elseif strcmp(type,'tri')

        prompt{end+1} = 'CSA Cf (pF)';
        def{end+1} = '0.2';

        prompt{end+1} = 'CSA I_reset (nA)';
        def{end+1} = '36';

    end
end

answ = inputdlg(prompt,'Parameter Values',[1 60],def);

% =============================
% PARSE PARAMETERS
% =============================
Cf = []; Rf = []; I_reset = 0;

tau_CR = zeros(1,n_CR);
R_cr   = zeros(1,n_CR);
C_cr   = zeros(1,n_CR);

tau_RC = zeros(1,n_RC);
R_rc   = zeros(1,n_RC);
C_rc   = zeros(1,n_RC);

idx = 1;

for i=1:length(labels)

    name = labels{i};
    type = types{i};

    % --- CSA ---
    if i==1

        if strcmp(type,'tau')

            tau_csa = str2double(answ{idx})*1e-9;
            idx = idx + 1;

            Cf = 0.2e-12;
            Rf = tau_csa / Cf;

        elseif strcmp(type,'rc')

            Rf = str2double(answ{idx})*1e3;
            Cf = str2double(answ{idx+1})*1e-12;
            idx = idx + 2;

            tau_csa = Rf * Cf;

        elseif strcmp(type,'tri')

            Cf = str2double(answ{idx})*1e-12;
            I_reset = str2double(answ{idx+1})*1e-9;
            idx = idx + 2;

            tau_csa = NaN;
            Rf = 0;
        end

    % --- CR STAGES ---
    elseif contains(name,'CR')

        k = sscanf(name,'CR%d');

        if strcmp(type,'tau')

            tau_CR(k) = str2double(answ{idx})*1e-9;
            idx = idx + 1;

            C_cr(k) = 0.2e-12;
            R_cr(k) = tau_CR(k)/C_cr(k);

        else

            R_cr(k) = str2double(answ{idx})*1e3;
            C_cr(k) = str2double(answ{idx+1})*1e-12;
            idx = idx + 2;

            tau_CR(k) = R_cr(k)*C_cr(k);
        end

    % --- RC STAGES ---
    elseif contains(name,'RC')

        k = sscanf(name,'RC%d');

        if strcmp(type,'tau')

            tau_RC(k) = str2double(answ{idx})*1e-9;
            idx = idx + 1;

            C_rc(k) = 0.2e-12;
            R_rc(k) = tau_RC(k)/C_rc(k);

        else

            R_rc(k) = str2double(answ{idx})*1e3;
            C_rc(k) = str2double(answ{idx+1})*1e-12;
            idx = idx + 2;

            tau_RC(k) = R_rc(k)*C_rc(k);
        end
    end
end


%% =============================
% BUILD PARAM STRUCT
%% =============================

params = struct();

% --- CSA ---
if csa_mode
    params.CSA.mode = 'exp';
    params.CSA.Cf   = Cf;
    params.CSA.Rf   = Rf;
else
    params.CSA.mode = 'tri';
    params.CSA.Cf   = Cf;
    params.CSA.I_reset = I_reset;
end

% --- FILTER CHAIN ---
params.shaper.chain = {};

% Add CR stages
for k = 1:n_CR
    stage.type = 'CR';
    stage.tau  = tau_CR(k);
    params.shaper.chain{end+1} = stage;
end

% Add RC stages
for k = 1:n_RC
    stage.type = 'RC';
    stage.tau  = tau_RC(k);
    params.shaper.chain{end+1} = stage;
end



%% =============================
% DISPLAY PARAMETERS
%% =============================

fprintf('\n=== ELECTRONICS PARAMETERS ===\n');

if strcmp(params.CSA.mode,'exp')
    fprintf('CSA: Exponential\n');
    fprintf('Cf = %.2e F, Rf = %.2e Ohm\n', Cf, Rf);
    fprintf('Tau_CSA = %.2f ns\n', tau_csa*1e9);
elseif strcmp(params.CSA.mode,'tri')
    fprintf('CSA: Triangular (constant current)\n');
    fprintf('Cf = %.2e F, I_reset = %.2e A\n', Cf, I_reset);
end

for k = 1:n_CR
    fprintf('CR%d: tau = %.2f ns (R=%.2e, C=%.2e)\n', ...
        k, tau_CR(k)*1e9, R_cr(k), C_cr(k));
end


for k = 1:n_RC
    fprintf('RC%d: tau = %.2f ns (R=%.2e, C=%.2e)\n', ...
        k, tau_RC(k)*1e9, R_rc(k), C_rc(k));
end

%% =============================
% TIME + INPUT
%% =============================

%data = load_preprocessed_series('HEXITEC_MHZ_CdTeFromAttemptedCZT_250micron_2mm_1000V_centre_both_%dns',35); %THIS NEEDS TO BE AN IMPORT LATER.
load("fastermaptest.mat");
dt = 1e-9 %CHANGE THIS TO A REQUEST FROM USER AT SOME POINT, OR IMPORT WITH MAP
params.dt = dt;

[Qmap_txyz, nx, ny, nz] = reshape_to_4D(data);

%% =============================
% TEST VOXEL
%% =============================

ix=round(nx/2); iy=round(ny/2); iz=round(nz/2);
Q = squeeze(Qmap_txyz(:,ix,iy,iz));
Q = Q / max(Q); % THIS LINE NEEDS REMOVING WHEN PRETTYMAN IS ALREADY SCALED BY CIE PROCESS, WHICH MAY BE SUBOPTIMAL (Q_max < 1).
Q_event = xray_to_charge(100, 4.4); %SCALE FOR NOW BASED ON A TEST PHOTON ENERGY OF 100 keV and work function of CdTe.
Q = Q * Q_event;

nt = length(Q);
t = (0:nt-1)*dt;

I = gradient(Q,dt);

%% =============================
% SIGNAL GENERATION (BOTH CSA MODES)
%% =============================

trace_length = 200; %Should be either a variable or defined  by user eventually
params.trace_length = trace_length * 1e-9; % seconds

if strcmp(params.CSA.mode,'exp')
V_exp = CSA_model(I,t,dt,params.CSA.Cf,params.CSA.Rf,'exp',0, trace_length);
elseif strcmp(params.CSA.mode,'tri')
V_tri = CSA_model(I,t,dt,params.CSA.Cf,0,'tri',params.CSA.I_reset, trace_length);
end
%{
V_exp = CSA_model(I,t,dt,Cf,Rf,'exp',0, trace_length);
V_tri = CSA_model(I,t,dt,Cf,0,'tri',I_reset, trace_length);
%}

if use_shaper
    if strcmp(params.CSA.mode,'exp')
    V_exp = apply_shaper(V_exp, params);
    end
    if strcmp(params.CSA.mode,'tri')
    V_tri = apply_shaper(V_tri, params);
    end
end


%{
if use_shaper
    V_exp = apply_shaper(V_exp,dt,tau_CR,n_RC,tau_RC);
    V_tri = apply_shaper(V_tri,dt,tau_CR,n_RC,tau_RC);
end
%}

%% =============================
% PLOT
%% =============================

figure(1);
hold off
if strcmp(params.CSA.mode,'exp')
plot([0:dt:trace_length/1e9-dt], V_exp*1000,'b','LineWidth',2); hold on;
legend('Exponential CSA');
end
if strcmp(params.CSA.mode,'tri')
plot([0:dt:trace_length/1e9-dt], V_tri*1000,'r--','LineWidth',2);
legend('Triangular CSA');
end
xlabel('Time (s)');
ylabel('Signal (mV)');
title('CSA + Shaper Comparison');
grid on;

%%
%{
figure(3);
hold off
plot([0:dt:trace_length/1e9-dt], exptau*1000,'b','LineWidth',2); hold on;
hold on
plot([0:dt:trace_length/1e9-dt], expRC*1000,'r--','LineWidth',2);
legend('exptau','expRC');
xlabel('Time (s)');
ylabel('Signal (mV)');
title('CSA + Shaper Comparison');
grid on;
%}


%%
load('VerificationData')
pairs = {
    tcsa_exp,  rcsa_exp;
    tcsa_tri, rcsa_tri;
    texp11,   rexp11;
    texp12,   rexp12;
    ttri11,   rtri11;
    ttri12,   rtri12
};

pair_names = {
    'tCSA exp vs RCSA exp';
    'TCSA tri vs RCSA tri';
    'Texp11 vs Rexp11';
    'Texp12 vs Rexp12';
    'Ttri11 vs Rtri11';
    'Ttri12 vs Rtri12'
};


for k = 1:size(pairs,1)
    figure(k);
    plot(pairs{k,1}, 'LineWidth', 1.5); hold on;
    plot(pairs{k,2}, '--','LineWidth', 1.5);

    title(pair_names{k});
    legend({'Signal 1','Signal 2'}, 'Location','best');
    grid on;
end



figure(7)
 plot(tcsa_exp, 'LineWidth', 1.5); hold on;
    plot(tcsa_tri, '--','LineWidth', 1.5);
