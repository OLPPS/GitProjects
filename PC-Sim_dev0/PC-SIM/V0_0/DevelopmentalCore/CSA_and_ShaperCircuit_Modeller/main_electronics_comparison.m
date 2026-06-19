clearvars -except exp_stored tri_stored exp_stored2 tri_stored2; clc; close all;

%% =============================
% CSA TYPE (CLICK SELECTION)
%% =============================

choice = questdlg('Select CSA Type:', ...
    'CSA Selection', ...
    'Exponential (RC)','Triangular (Constant Current)', ...
    'Exponential (RC)');

if strcmp(choice,'Exponential (RC)')
    csa_mode = 1;
else
    csa_mode = 2;
end

%% =============================
% CSA PARAMETERS
%% =============================

if csa_mode == 1

    choice = questdlg('Define CSA using:', ...
        'CSA Input Mode', ...
        'Tau only','R & C','Tau only');

    if strcmp(choice,'Tau only')

        answ = inputdlg({'Enter tau_CSA (ns)'}, ...
            'CSA Tau',[1 60],{'10000'});

        tau_csa = str2double(answ{1}) * 1e-9;

        % Assign default capacitance
        Cf = 0.2e-12;   % default
        Rf = tau_csa / Cf;

    else

        answ = inputdlg({'Enter Cf (pF)','Enter Rf (MOhm)'}, ...
            'CSA R&C',[1 60],{'0.2','50'});

        Cf = str2double(answ{1}) * 1e-12;
        Rf = str2double(answ{2}) * 1e6;

        tau_csa = Rf * Cf;
    end

    I_reset = 0;

else

    answ = inputdlg({'Enter Cf (pF)','Enter reset current (nA)'}, ...
        'CSA Triangular',[1 60],{'0.2','36'});

    Cf      = str2double(answ{1}) * 1e-12;
    I_reset = str2double(answ{2}) * 1e-9;

    Rf = NaN;
    tau_csa = NaN;

end

%% =============================
% SHAPER ON/OFF
%% =============================

choice = questdlg('Apply shaping?', ...
    'Shaper Selection', ...
    'Yes','No','Yes');

use_shaper = strcmp(choice,'Yes');

%% =============================
% CR STAGE
%% =============================

use_CR = 0;
tau_CR = 0;
R_cr = NaN; C_cr = NaN;

if use_shaper

    choice = questdlg('Use CR stage?', ...
        'CR Selection', ...
        'Yes','No','Yes');

    use_CR = strcmp(choice,'Yes');

    if use_CR

        choice = questdlg('Define CR using:', ...
            'CR Input Mode', ...
            'Tau','R & C','Tau');

        if strcmp(choice,'Tau')
            answ = inputdlg({'Enter tau_CR (ns)'},'CR Tau',[1 60],{'20'});
            tau_CR = str2double(answ{1})*1e-9;

            % equivalent R/C (choose default C)
            C_cr = 0.2e-12;
            R_cr = tau_CR / C_cr;

        else
            answ = inputdlg({'Enter R (kOhm)','Enter C (pF)'}, ...
                'CR R&C',[1 60],{'0.2','100'});

            R_cr = str2double(answ{1})*1e3;
            C_cr = str2double(answ{2})*1e-12;
            tau_CR = R_cr * C_cr;
        end
    end
end

%% =============================
% RC STAGES
%% =============================

n_RC = 0;

n_RC = 0;
tau_RC = [];
R_rc = []; 
C_rc = [];

if use_shaper

    answ = inputdlg({'Number of RC stages (n)'}, ...
        'RC Stages',[1 60],{'1'});

    n_RC = str2double(answ{1});

    tau_RC = zeros(1,n_RC);
    R_rc   = zeros(1,n_RC);
    C_rc   = zeros(1,n_RC);

    if n_RC > 0

        % --- choose input mode once ---
        choice = questdlg('Define RC stages using:', ...
            'RC Input Mode','Tau','R & C','Tau');

        % --- build prompts dynamically ---
        prompt = cell(1, n_RC * (strcmp(choice,'Tau') + 2*strcmp(choice,'R & C')));
        definput = cell(size(prompt));

        idx = 1;

        for k = 1:n_RC

            if strcmp(choice,'Tau')

                prompt{idx} = sprintf('Stage %d: tau (ns)', k);
                definput{idx} = '20';
                idx = idx + 1;

            else

                prompt{idx} = sprintf('Stage %d: R (kOhm)', k);
                definput{idx} = '100';
                idx = idx + 1;

                prompt{idx} = sprintf('Stage %d: C (pF)', k);
                definput{idx} = '0.2';
                idx = idx + 1;

            end
        end

        % --- single input dialog ---
        answ = inputdlg(prompt, 'RC Stage Parameters', [1 60], definput);

        % --- parse input ---
        idx = 1;

        for k = 1:n_RC

            if strcmp(choice,'Tau')

                tau_RC(k) = str2double(answ{idx}) * 1e-9;

                % assign equivalent R,C
                C_rc(k) = 0.2e-12;   % default
                R_rc(k) = tau_RC(k) / C_rc(k);

                idx = idx + 1;

            else

                R_rc(k) = str2double(answ{idx}) * 1e3;
                C_rc(k) = str2double(answ{idx+1}) * 1e-12;

                tau_RC(k) = R_rc(k) * C_rc(k);

                idx = idx + 2;
            end
        end

    end
end


%% =============================
% DISPLAY PARAMETERS
%% =============================

fprintf('\n=== ELECTRONICS PARAMETERS ===\n');

if csa_mode == 1
    fprintf('CSA: Exponential\n');
    fprintf('Cf = %.2e F, Rf = %.2e Ohm\n', Cf, Rf);
    fprintf('Tau_CSA = %.2f ns\n', tau_csa*1e9);
else
    fprintf('CSA: Triangular (constant current)\n');
    fprintf('Cf = %.2e F, I_reset = %.2e A\n', Cf, I_reset);
end

if use_CR
    fprintf('CR: tau = %.2f ns (R=%.2e, C=%.2e)\n', ...
        tau_CR*1e9, R_cr, C_cr);
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

trace_length = 200;

V_exp = CSA_model(I,t,dt,Cf,Rf,'exp',0, trace_length);
V_tri = CSA_model(I,t,dt,Cf,0,'tri',I_reset, trace_length);

if use_shaper
    V_exp = apply_shaper(V_exp,dt,use_CR,tau_CR,n_RC,tau_RC);
    V_tri = apply_shaper(V_tri,dt,use_CR,tau_CR,n_RC,tau_RC);
end

%% =============================
% PLOT
%% =============================

figure(1);
hold off
plot([0:dt:trace_length/1e9-dt], V_exp*1000,'b','LineWidth',2); hold on;
plot([0:dt:trace_length/1e9-dt], V_tri*1000,'r--','LineWidth',2);

legend('Exponential CSA','Triangular CSA');
xlabel('Time (s)');
ylabel('Signal (mV)');
title('CSA + Shaper Comparison');
grid on;

%{
clear; clc; close all;

%% =============================
% GUI INPUT
%% =============================

% --- CSA TYPE ---

hoice = questdlg('Select CSA Type:', ...
                  'CSA Selection', ...
                  'Exponential (RC)','Triangular','Exponential (RC)');

if strcmp(choice,'Exponential (RC)')
    csa_mode = 1;
else
    csa_mode = 2;
end


% --- CSA PARAMETERS ---
if csa_mode == 1
    prompt = {'Enter Cf (pF)','Enter Rf (MOhm)'};
    def = {'0.2','50'};
    answ = inputdlg(prompt,'CSA Exponential',[1 60],def);
    Cf = str2double(answ{1})*1e-12;
    Rf = str2double(answ{2})*1e6;
    I_reset = 0;

elseif csa_mode == 2
    prompt = {'Enter Cf (pF)','Enter reset current (uA)'};
    def = {'0.2','5'};
    answ = inputdlg(prompt,'CSA Triangular',[1 60],def);
    Cf = str2double(answ{1})*1e-12;
    I_reset = str2double(answ{2})*1e-6;
    Rf = 0;
end

% --- SHAPER ON/OFF ---
choice = questdlg('Apply shaping?', ...
                  'Shaper Selection', ...
                  'Yes','No','Yes');

if strcmp(choice,'Yes')
    use_shaper = 1;
else
    use_shaper = 0;
end


use_CR = 0;
n_RC = 0;

if use_shaper == 1

    % --- CR ---
    
choice = questdlg('Use CR stage?', ...
                  'CR Selection', ...
                  'Yes','No','Yes');

use_CR = strcmp(choice,'Yes');

    if use_CR
        prompt = {'Enter tau_CR (ns)'};
        answ = inputdlg(prompt,'CR tau',[1 60],{'20'});
        tau_CR = str2double(answ{1})*1e-9;
    else
        tau_CR = 0;
    end

    % --- RC ---
    prompt = {'Number of RC stages (n)'};
    answ = inputdlg(prompt,'RC stages',[1 60],{'1'});
    n_RC = str2double(answ{1});

    tau_RC = zeros(1,n_RC);

    for k = 1:n_RC
        prompt = {sprintf('tau_RC%d (ns)',k)};
        answ = inputdlg(prompt,'RC tau',[1 60],{'20'});
        tau_RC(k) = str2double(answ{1})*1e-9;
    end
else
    tau_CR = 0;
    tau_RC = [];
end

%% =============================
% TIME + INPUT SIGNAL
%% =============================

dt = 1e-9;
t  = 0:dt:200e-9;
nt = length(t);

% TEST INPUT FROM DETECTOR
Q = exp(-(t-20e-9).^2/(2*(5e-9)^2));
Q = Q / max(Q);

I = gradient(Q,dt);

%% =============================
% RUN BOTH CSA MODES FOR COMPARISON
%% =============================

V_exp = CSA_model(I,t,dt,Cf,50e6,'exp',0); % reference exp
V_tri = CSA_model(I,t,dt,Cf,0,'tri',I_reset);

if use_shaper

    V_exp = apply_shaper(V_exp,dt,use_CR,tau_CR,n_RC,tau_RC);
    V_tri = apply_shaper(V_tri,dt,use_CR,tau_CR,n_RC,tau_RC);

end

%% =============================
% PLOT
%% =============================

figure;

plot(t, V_exp,'b','LineWidth',2); hold on;
plot(t, V_tri,'r--','LineWidth',2);

legend('Exponential CSA','Triangular CSA');
xlabel('Time (s)');
ylabel('Signal');
title('CSA Comparison (Exponential vs Triangular)');
grid on;

%}