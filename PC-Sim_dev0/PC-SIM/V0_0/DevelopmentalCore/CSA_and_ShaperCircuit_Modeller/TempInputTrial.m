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
    % ===== MIXED MODE (your original version, slightly cleaned) =====

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
% PARSE PARAMETERS (UNCHANGED)
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

            Rf = str2double(answ{idx})*1e6;
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