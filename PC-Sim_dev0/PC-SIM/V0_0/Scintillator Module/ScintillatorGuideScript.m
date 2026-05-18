%This script is for command line based input
clear
clc

%Conversion properties for excited site to optical photons - Light Yield
EnergyDependance = input("Is there energy dependence? (Y/N): ","s");

%%%%%%%%CheckIfValidInput

switch EnergyDependance
    case "N"
disp('Variability in number of photons produced is modelled as a normal distribution. Please input: ');
Exci2Opti_mean = input("Average Light Yield Per keV: ");
Exci2Opti_std = input('Standard deviation of the distribution: ');
    case "Y"
%TO DO
end

%Dispersion of emissions in time - excited site decay
disp("The decay of optical photons is modelled with exponential decays.");
NumOfExps = input("How many independent decays would you like to use?");
%%%%%%%%CheckIfValidInput
ExponentialDecayConstants_ = zeros(NumOfExps,1);

for ModelledDecayNum = 1:NumOfExps
fprintf("For exponential decay number %i", ModelledDecayNum)
    ExponentialDecayConstants(ModelledDecayNum) = input(' what is the decay constant (in ns)?: ');
end

% Fraction of Light reaching SiPM - Light survival
%%%%%THIS COULD BE CHANGED IN FUTURE TO BE SPATIALLY DEPENDENT
disp("Some fraction of optical photons will not reach the SiPM due to crosstalk, escape or capture by septa.");
disp("The surviving fraction of photons should thus be between 0 (complete loss) and 1 (no losses).");
disp("Note that fill factor of SPADs in the SiPM is NOT to be included in this step, just lossess associated with pre-SiPM transmission.")
SurvivalFraction = input("For this model, what fraction of the emitted optical photons are expected to reach the SiPM?: ");
%%%%%%%%CheckIfValidInput


%Generate Scintillator parameters file here - ScintCrystalParams
save("ScintCrystalParams","EnergyDependance","Exci2Opti_mean","Exci2Opti_std","ExponentialDecayConstants","SurvivalFraction")



%%
%This script is for GUI based input

%{
% --- Energy dependence (Yes/No) ---
choice = questdlg('Is there energy dependence?', ...
                  'Energy Dependence', ...
                  'Yes', 'No', 'No');
EnergyDependance = choice;   % Returns 'Yes' or 'No'

% --- Numeric inputs ---
prompt = { ...
    'Average Number Of Photons Produced Per Deposition:', ...
    'Statistical Spread Of Average (Mean):', ...
    'Statistical Spread Of Average (Standard deviation):', ...
    'Statistical Spread Of Average (Mean):'};

title = 'Input Parameters';
dims = [1 50];

answer = inputdlg(prompt, title, dims);

AvgPhotConvRatio = str2double(answer{1});
Exci2Opti_mean = str2double(answer{2});
Exci2Opti_std = str2double(answer{3});

E = str2double(answer{4});


% Create a Yes/No dialog box
choice = questdlg('Do you want to continue?', ... % Question
                  'Confirmation', ...             % Dialog title
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