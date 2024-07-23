% ChemControl_cbm_eval_param

% This is an interactive script---execute it step-by-step.
% It fits a serious of computational reinforcement learning models using the CBM toolbox.
% and evaluates them.
%
% EXPLANATION OF SETTINGS:
% dirs.root               = string, root directory of project.
%
% OUTPUTS:
% no outputs, just plots.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.
% Inspired by code from Algermissen, J. et al. 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel
clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Set directories:

% Directories:
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.results = fullfile(dirs.root, 'Log', 'Behavior', 'Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
dirs.models         = fullfile(dirs.root, 'models');

% Add paths:
addpath('/home/control/samnel/Documents/MATLAB/cbm-master/codes'); % CBM toolbox
addpath(fullfile(dirs.root, 'behavioral_study', 'scripts', 'matlab_scripts', 'RescorlaWagnerModel', 'models')); % models

% Set transforms and parameter names:
param_names{1} = {'\epsilon', '\rho'}; % Model
transform{1} = {'sigmoid', 'exp'};

param_names{2} = {'\epsilon', '\rho', 'goBias'}; % GoBiasModel
transform{2} = {'sigmoid', 'exp', '@(x) x'};

param_names{3} = {'\epsilon', '\rho', 'goBias', '\pi'}; % FixedPavlovModel
transform{3} = {'sigmoid', 'exp', '@(x) x', '@(x) x'};

param_names{4} = {'\epsilon', '\rho', 'goBias', '\pi'}; % DynamicPavlovModel
transform{4} = {'sigmoid', 'exp', '@(x) x', '@(x) x'};

param_names{5} = {'\epsilon', '\rho', 'goBias', '\omega'}; % FixedOmegaModel
transform{5} = {'sigmoid', 'exp', '@(x) x', 'sigmoid'};

param_names{6} = {'\epsilon', '\rho', 'goBias', '\alpha', '\kappa', '\slope'}; % DynamicOmega1Model
transform{6} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid', 'exp'};

param_names{7} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{7} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

param_names{8} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{8} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

param_names{9} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', 'lapse bounds'}; % DynamicOmega2Model
transform{9} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid', 'sigmoid'};

param_names{10} = {'\epsilon', '\rho', 'goBias','\beta_{\Omega}', '\thres_{\Omega}', 'reward info'}; % DynamicOmega2Model
transform{10} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid', '@(x) x'};

param_names{11} = {'\epsilon', '\rho', 'goBias'}; % Model 11
transform{11} = {'sigmoid', 'exp', '@(x) x'};

param_names{12} = {'\epsilon', '\rho', 'goBias', '\alpha_{up}', '\alpha_{down}'}; % Model 12
transform{12} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid'};

param_names{13} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{13} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid', 'sigmoid'};

param_names{14} = {'\epsilon', '\beta_{softmax}', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{14} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid', 'sigmoid'};

param_names{15} = {'\epsilon', '\beta_{softmax}', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{15} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid', 'sigmoid'};

param_names{16} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{16} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid', 'sigmoid'};

param_names{17} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{17} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

param_names{18} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{18} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

param_names{19} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{19} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

param_names{20} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{20} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

param_names{21} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{21} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', 'scaledSigmoid'};

% ----------------------------------------------------------------------- %
%% 00b) Select model:

% Get list of .m files in dir.models
modelFiles = dir(fullfile(dirs.models, '*.m'));

% Set nMod to the number of .m files
nMod = length(modelFiles);

selMod = 7;

fprintf("Selected model is no. %02d\n", selMod)

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 01) LOAD AND SAVE ALL DATA:
% 01a) Alternative A: Load model fitted with LAP:

parType = 'lap';
dataType = 'allData';
fname_mod = cell(nMod, 1);

for iMod = 1:nMod
    fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, dataType));
end

fprintf("Load model %d fit with LAP\n", selMod);
fname       = load(fname_mod{selMod});
cbm         = fname.cbm;
subParam    = cbm.output.parameters;
nSub        = size(subParam, 1);
nParam      = size(subParam, 2);

% ----------------------------------------------------------------------- %
% 01b) Alternative B: Load model fitted with HBI:

parType = 'hbi';

% load
fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d_%s.mat', selMod, dataType));
fprintf("Load model %d fit with HBI\n", selMod);
load(fname_hbi)

% extract
fprintf("Extract model %d fit with HBI\n", selMod);
groupParam      = cbm.output.group_mean{:};
subParam        = cbm.output.parameters{:};
nSub            = size(subParam, 1);
nParam          = size(subParam, 2);



% ----------------------------------------------------------------------- %
% 01c) Transform group and subject parameters appropriately:
% Define sigmoid function
sigmoid = @(x) 1 ./ (1 + exp(-x));

fprintf('Transform parameters of model %d\n', selMod);

for iParam = 1:nParam
    % Get the transformation function for the current parameter
    transformFunc = transform{selMod}{iParam};
    transformFunc = str2func(transformFunc);
    groupParam(iParam) = transformFunc(groupParam(iParam));
    subParam(:, iParam) = transformFunc(subParam(:, iParam));
end

% ----------------------------------------------------------------------- %
% 01d) Save (either LAP or HBI):

fprintf('Save transformed group and subject level parameters under %s, model %02d\n', ...
    parType, selMod);

fileName = fullfile(eval(sprintf('dirs.%s', parType)), sprintf('CBM_%s_M%02d_%s_groupParameters.csv', ...
    parType, selMod, dataType));

csvwrite(fileName, groupParam);

fileName = fullfile(eval(sprintf('dirs.%s', parType)), sprintf('CBM_%s_M%02d_%s_subjectParameters.csv', ...
    parType, selMod, dataType));

csvwrite(fileName, subParam);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% 01) GROUP-LEVEL PARAMETERS.
% Mean per proup level parameter:

fprintf('\nGroup-level parameters for model %d:\n', selMod);

for iParam = 1:size(groupParam, 2)
    x = groupParam(:, iParam);
    fprintf('Parameter %d: M = %.02f\n', iParam, mean(x));
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% 01) SUBJECT-LEVEL PARAMETERS.
% Mean, SD, range per subject-level parameter:

fprintf('\nSubject-level parameters for model %d\n', selMod);
for iParam = 1:size(subParam, 2)
    x = subParam(:, iParam);
    fprintf('Parameter %d: M = %.03f, SD = %.03f, range %.03f - %.03f \n', ...
        iParam, mean(x), std(x), min(x), max(x));
end

% ----------------------------------------------------------------------- %
% 01) Find extreme values:

fprintf('Find subjects with extreme values for model %d\n', selMod);
for iParam = 1:size(subParam, 2)
    x       = subParam(:, iParam);
    xmin    = min(x);
    xminidx = find(x == xmin);
    xmax    = max(x);
    xmaxidx = find(x == xmax);
    fprintf('Parameter %d: min = %.05f for subjects %s, max = %.02f for subjects %s \n', ...
        iParam, xmin, num2str(xminidx), xmax, num2str(xmaxidx));
end

% ----------------------------------------------------------------------- %
% 01) Sign of parameter:

fprintf('Determine percentage negative parameters for model %d\n', selMod);
for iParam = 1:size(subParam, 2)
    x       = subParam(:, iParam);
    xNeg    = sum(x < 0); % how many negative
    xIdx    = find(x < 0)'; % who negative
    fprintf('Parameter %d: negative for %d subjects: %s \n', ...
        iParam, xNeg, num2str(xIdx));
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% 01) PLOTS.
% Plot parameters with cbm function:
% Parameter names:

% Model names:
model_names = {'M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21'};
% Create output name:
modVec = 1:nMod;
fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s_%s.mat', num2str(modVec, '_%02d'), dataType));
d = load(fname_hbi);
c = d.cbm;
freq = c.output.model_frequency;
[~, k] = max(freq);

%cbm_hbi_plot(fname_hbi, model_names, param_names{k}, transform{k});

%%
control_types = {'allData', 'hc', 'lc'};
parTypes = {'hbi'};


for ctype = control_types
    ctype = ctype{1}; % Extract string from cell
    
    for parType = parTypes
        parType = parType{1}; % Extract string from cell

        % Load model fitted with LAP or HBI
        if strcmp(parType, 'lap')
            fname_mod = cell(nMod, 1);
            for iMod = 1:nMod
                fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, ctype));
            end

            fprintf("Load model %d fit with LAP (%s)\n", selMod, ctype);
            fname = load(fname_mod{selMod});
            cbm = fname.cbm;
            subParam = cbm.output.parameters;
            nSub = size(subParam, 1);
            nParam = size(subParam, 2);
        else
            fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod_%02d_%s.mat', selMod, ctype));
            fprintf("Load model %d fit with HBI (%s)\n", selMod, ctype);
            load(fname_hbi)

            fprintf("Extract model %d fit with HBI (%s)\n", selMod, ctype);
            groupParam = cbm.output.group_mean{:};
            subParam = cbm.output.parameters{:};
            nSub = size(subParam, 1);
            nParam = size(subParam, 2);
        end

        % Transform group and subject parameters appropriately
        sigmoid = @(x) 1 ./ (1 + exp(-x));

        fprintf('Transform parameters of model %d (%s)\n', selMod, ctype);
        for iParam = 1:nParam
            transformFunc = transform{selMod}{iParam};
            transformFunc = str2func(transformFunc);
            
            groupParam(iParam) = transformFunc(groupParam(iParam));
            subParam(:, iParam) = transformFunc(subParam(:, iParam));
        end

        % Save transformed group and subject level parameters
        fprintf('Save transformed group and subject level parameters under %s, model %02d (%s)\n', ...
            parType, selMod, ctype);

        fileName = fullfile(eval(sprintf('dirs.%s', parType)), sprintf('CBM_%s_M%02d_%s_groupParameters.csv', ...
            parType, selMod, ctype));
        csvwrite(fileName, groupParam);

        fileName = fullfile(eval(sprintf('dirs.%s', parType)), sprintf('CBM_%s_M%02d_%s_subjectParameters.csv', ...
            parType, selMod, ctype));
        csvwrite(fileName, subParam);

        % Group-level parameters
        fprintf('\nGroup-level parameters for model %d (%s):\n', selMod, ctype);
        for iParam = 1:size(groupParam, 2)
            x = groupParam(:, iParam);
            fprintf('Parameter %d: M = %.02f\n', iParam, mean(x));
        end

        % Subject-level parameters
        fprintf('\nSubject-level parameters for model %d (%s):\n', selMod, ctype);
        for iParam = 1:size(subParam, 2)
            x = subParam(:, iParam);
            fprintf('Parameter %d: M = %.03f, SD = %.03f, range %.03f - %.03f \n', ...
                iParam, mean(x), std(x), min(x), max(x));
        end

        % Find extreme values
        fprintf('Find subjects with extreme values for model %d (%s):\n', selMod, ctype);
        for iParam = 1:size(subParam, 2)
            x = subParam(:, iParam);
            xmin = min(x);
            xminidx = find(x == xmin);
            xmax = max(x);
            xmaxidx = find(x == xmax);
            fprintf('Parameter %d: min = %.05f for subjects %s, max = %.02f for subjects %s \n', ...
                iParam, xmin, num2str(xminidx), xmax, num2str(xmaxidx));
        end

        % Sign of parameter
        fprintf('Determine percentage negative parameters for model %d (%s):\n', selMod, ctype);
        for iParam = 1:size(subParam, 2)
            x = subParam(:, iParam);
            xNeg = sum(x < 0); % how many negative
            xIdx = find(x < 0)'; % who negative
            fprintf('Parameter %d: negative for %d subjects: %s \n', ...
                iParam, xNeg, num2str(xIdx));
        end

        % Plots
        % Create output name
        fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s_%s.mat', num2str(modVec, '_%02d'), ctype));
        d = load(fname_hbi);
        c = d.cbm;
        freq = c.output.model_frequency;
        [~, k] = max(freq);

        cbm_hbi_plot(fname_hbi, model_names, param_names{k}, transform{k});
    end
end

% END OF FILE.