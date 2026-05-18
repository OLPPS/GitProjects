%This script is for command line based input
clear
clc

%Photon detection efficiency - PDE
PDE = input("Of the photons which make it to the SiPM, what fraction are detected (what is the PDE, from 0 - 1)?: ");
%%%%%%PDE NEEDS TO EVENTUALLY HAVE A FUNCTION OF LAMBDA OPTION
disp("SiPM may be square or rectangular. Which direction is X and which is Y does not matter for this step however.");
SPADsInX = input("How many SPADs are in the x-direction of the SPAD?: ");
SPADsInY = input("How many SPADs are in the y-direction of the SPAD?: ");
%TotalNumOfSPADs = SPADsInX * SPADsInY;


%Dispersion of recharges in time is modelled as
disp("The recharge of a given SPAD is modelled with exponential decays.");
NumOfExps = input("How many independent decays would you like to use in the model?");
%%%%%%%%CheckIfValidInput
ExponentialDecayConstants_SPADsRecharge = zeros(NumOfExps,1);

for ModelledDecayNum = 1:NumOfExps
    fprintf("For exponential decay number %i", ModelledDecayNum);
    ExponentialDecayConstants_SPADsRecharge(ModelledDecayNum) = input(' what is the decay constant (in ns)?: ');
end

%Cross-talk parameters
FractionalCrossTalk = input("What fraction of the SPAD triggers would be expected to trigger cross-talk in an another SPAD (in the same SiPM)? (0 - 1): ");

%Should SPADs vary in their gain
disp("This tool allows the gain of SPADs within an SiPM to vary.")
disp("These variations are modelled with a normal distribution.")
SPADsVariations = input("Do you want to model varations in SPAD gain? (Y/N)","s");
%%%%%%%%CheckIfValidInput
%%%%% IF YES, USE NOISE MODELLER LATER
SPADsGainNoise_Mean = 1;
SPADsGainNoise_STD = 0;
SPADsGainNoise_Cap = 0;
CapVariation = "N";
if SPADsVariations == "Y"
    SPADsGainNoise_Mean = input('What is the mean gain the SPADs should have (default is 1): ');
    SPADsGainNoise_STD = input('What is the standard deviation in gain the SPADs should have (default is 1): ');
    disp("Do you want to cap the range of gains that the SPADs can take? (Y/N)")
    CapVariation = input("Note that gains will not go below 0 by default.: ");
    if CapVariation == "Y"
        SPADsGainNoise_Cap = input("What is the maximum gain variation allowed (in sigma)?: ");
    end
end

%Generate SiPM parameters file here - SiPMParams
save("SiPMParams","PDE","SPADsInX","SPADsInY","ExponentialDecayConstants_SPADsRecharge","FractionalCrossTalk","SPADsVariations","SPADsGainNoise_Mean", "SPADsGainNoise_STD", "SPADsGainNoise_Cap", "CapVariation");

%Calculate








%%
%This script is for GUI based input
clear
clc

% --- Recall code or new SiPM? --- 
RecallSiPM = questdlg('Do you have an SiPM recall code?', 'SiPM recall', 'Yes', 'No', 'No');

if RecallSiPM == 'Yes'
prompt = { ...
    'Please enter SiPM Recall Code:'}

TITLE = 'Physical SiPM parameters';
dims = [1 100];

encodedparams = ForcedFilledDialogueBox(prompt, TITLE, dims);
    
else
% --- Is response wavelength independent? ---
ResponseLinearity = questdlg('Do you want to vary parameters as a function of optical wavelength?', 'Response linearity', 'Yes', 'No', 'No');

    
if RecallSiPM == 'No' & ResponseLinearity == 'No'
% --- Physical SiPM parameters ---

   
prompt = { ...
    'Photon detection efficiency (what FRACTION of optical photons reaching a SPAD will cause it to fire):', ...
    'SPADs in X direction:', ...
    'SPADs in Y direction:' ...
    'Geometric filling (what FRACTION of the SiPM cross-section is SPADs as opposed to spacers):', ...
    'SPAD crosstalk probabilty (as a fraction):'}

TITLE = 'Physical SiPM parameters';
dims = [1 50];

answer = ForcedFilledDialogueBox(prompt, TITLE, dims);

PDE = str2double(answer{1});
SPADsInX = str2double(answer{2});
SPADsInY = str2double(answer{3});
FillFactor = str2double(answer{4});
FractionalCrossTalk = str2double(answer{5});

end



prompt = {"The recharge of a given SPAD is modelled with exponential decays. How many independent decays would you like to use in the model?"};

TITLE = ['SPAD discharge model parameters'];
dims = [1 50];

answer = inputdlg(prompt, TITLE, dims);

NumOfExps = str2double(answer{1});



% Preallocate prompt cell array
prompt = cell(NumOfExps, 1);

% Generate variable number of prompts
for ii = 1:NumOfExps
    prompt{ii} = sprintf('Decay constant %d is (in ns):', ii);
end

TITLE = 'Input decay constants';
dims = [1 50];

answer = inputdlg(prompt, TITLE, dims);

%Preallocate decya constants array
ExponentialDecayConstants_SPADsRecharge = zeros(NumOfExps,1);

%Fill array
for ii = 1:NumOfExps
ExponentialDecayConstants_SPADsRecharge(ii) = str2double(answer{ii});
end


% --- Recall code or new SiPM? --- 
SPADsVariations = questdlg('This tool allows the gain of SPADs within an SiPM to vary. These variations are modelled with a normal distribution. Do you want to model varations in SPAD gain?', 'SPAD gain variations', 'Yes', 'No', 'No');

SPADsGainNoise_Mean = 1;
SPADsGainNoise_STD = 0;
SPADsGainNoise_Cap = 0;
CapVariation = "N";
if SPADsVariations == "Yes"
    prompt = {'What is the mean gain the SPADs should have (default is 1):', ...
        'What is the standard deviation in gain the SPADs should have (default is 1):' ...
        }
        
TITLE = 'Input SPAD variation parameters';
dims = [1 50];

answer = inputdlg(prompt, TITLE, dims);

    SPADsGainNoise_Mean = str2double(answer{1});
    SPADsGainNoise_STD  = str2double(answer{2});

    CapVariation = questdlg("Do you want to cap the range of gains that the SPADs can take? Note that gains will not go below 0 by default.",'SPAD gain variations cap', 'Yes', 'No', 'No');
    if CapVariation == "Yes"
        prompt = {"What is the maximum gain variation allowed (in sigma)?"}
        TITLE = 'Set SPAD gain variation cap';
        sims = [1 50];

        SPADsGainNoise_Cap = inputdlg(prompt, TITLE, dims);
    end
end



%Calculate SiPM recall code

TAULIST = '';
for N = 1:NumOfExps
TAULIST = strcat(TAULIST,num2str(ExponentialDecayConstants_SPADsRecharge(N)),'n');
end

switch SPADsVariations
    case 'No'
SPADSVARASNUM = '0';
    case 'Yes'
SPADSVARASNUM = '1';
end


switch CapVariation
    case 'No'
CAPVARASNUM = '0';
    case 'Yes'
CAPVARASNUM = '1';
end


RecallCode_SiPM = strcat('PDE', num2str(PDE), 'XSP', num2str(SPADsInX), 'YSP', num2str(SPADsInY), 'DTAU1', TAULIST, 'CTF', num2str(FractionalCrossTalk), 'VARSP', SPADSVARASNUM, 'MNSP', num2str(SPADsGainNoise_Mean), 'STDSP', num2str(SPADsGainNoise_STD), 'CAPSP', num2str(SPADsGainNoise_Cap), 'VARCAP', CAPVARASNUM) ;

%Generate SiPM parameters file here - SiPMParams
encodedparams = decode_RecallCode_SiPM(RecallCode_SiPM);
end
save(strcat("SiPMParams_",RecallCode_SiPM),encodedparams);



%%
% Create a Yes/No dialog box
choice = questdlg('Do you want to continue?', ... % Question
                  'Confirmation', ...             % Dialog TITLE
                  'Yes', 'No', 'No');             % Button labels & default

% Handle the response
switch choice
    case 'Yes'
        disp('User chose YES.');
    case 'No'
        disp('User chose NO.');
    otherwise
        disp('User closed the dialog.');
end
%}

function answer = ForcedFilledDialogueBox(prompt, TITLE, dims)
while true
answer = inputdlg(prompt, TITLE, dims);



   % User pressed Cancel
    if isempty(answer)
        error('User cancelled input.');
    end

    % Check for empty entries
    if all(~cellfun(@isempty, answer))
        break
    else
        uiwait(warndlg('Please fill in all values before continuing.', ...
                       'Missing input'));
    end
end
end


%function (PDE,SPADsInX,SPADsInY,ExponentialDecayConstants_SPADsRecharge,FractionalCrossTalk,SPADsVariations,SPADsGainNoise_Mean, SPADsGainNoise_STD, SPADsGainNoise_Cap, CapVariation) = decodeRecallSCode_SiPM
function params = decode_RecallCode_SiPM(RecallCode_SiPM)
% Decode RecallCode_SiPM string into parameter values
%
% Output:
%   params : struct containing all decoded parameters

% ---------------------------
% 1. Extract scalar parameters
% ---------------------------
tokens = regexp(RecallCode_SiPM, ...
    ['PDE(?<PDE>[\d\.]+)' ...
     'XSP(?<XSP>\d+)' ...
     'YSP(?<YSP>\d+)' ...
     'DTAU1(?<TAU>[\dn\.]+)' ...
     'CTF(?<CTF>[\d\.]+)' ...
     'VARSP(?<VARSP>[01])' ...
     'MNSP(?<MNSP>[\d\.]+)' ...
     'STDSP(?<STDSP>[\d\.]+)' ...
     'CAPSP(?<CAPSP>[\d\.]+)' ...
     'VARCAP(?<VARCAP>[\d\.]+)'], ...
     'names');

if isempty(tokens)
    error('RecallCode_SiPM format not recognized.');
end

% ---------------------------
% 2. Convert to numeric values
% ---------------------------
params.PDE                      = str2double(tokens.PDE);
params.SPADsInX                 = str2double(tokens.XSP);
params.SPADsInY                 = str2double(tokens.YSP);
params.FractionalCrossTalk      = str2double(tokens.CTF);
params.SPADsGainNoise_Mean      = str2double(tokens.MNSP);
params.SPADsGainNoise_STD       = str2double(tokens.STDSP);
params.SPADsGainNoise_Cap       = str2double(tokens.CAPSP);
params.CapVariation             = str2double(tokens.VARCAP);

% ---------------------------
% 3. Decode SPAD variation flag
% ---------------------------
if strcmp(tokens.VARSP,'1')
    params.SPADsVariations = 'Yes';
else
    params.SPADsVariations = 'No';
end

% ---------------------------
% 4. Decode TAU list
% ---------------------------
% Original format: e.g. "10n25n50n"
tauStrings = regexp(tokens.TAU, '([\d\.]+)n', 'tokens');
params.ExponentialDecayConstants_SPADsRecharge = ...
    cellfun(@(x) str2double(x{1}), tauStrings);

end