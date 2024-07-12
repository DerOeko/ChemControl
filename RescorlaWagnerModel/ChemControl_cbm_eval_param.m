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
transform{1} = {@sigmoid, @exp};

param_names{2} = {'\epsilon', '\rho', 'goBias'}; % GoBiasModel
transform{2} = {@sigmoid, @exp, @(x) x};

param_names{3} = {'\epsilon', '\rho', 'goBias', '\pi'}; % FixedPavlovModel
transform{3} = {@sigmoid, @exp, @(x) x, @(x) x};

param_names{4} = {'\epsilon', '\rho', 'goBias', '\pi'}; % DynamicPavlovModel
transform{4} = {@sigmoid, @exp, @(x) x, @(x) x};

param_names{5} = {'\epsilon', '\rho', 'goBias', '\omega'}; % FixedOmegaModel
transform{5} = {@sigmoid, @exp, @(x) x, @sigmoid};

param_names{6} = {'\epsilon', '\rho', 'goBias', '\alpha', '\kappa', '\slope'}; % DynamicOmega1Model
transform{6} = {@sigmoid, @exp, @(x) x, @sigmoid, @sigmoid, @exp};

param_names{7} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{7} = {@sigmoid, @exp, @(x) x, @sigmoid, @exp, @scaledSigmoid};

param_names{8} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'}; % DynamicOmega2Model
transform{8} = {@sigmoid, @exp, @(x) x, @sigmoid, @exp, @scaledSigmoid};

param_names{9} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', 'lapse bounds'}; % DynamicOmega2Model
transform{9} = {@sigmoid, @exp, @(x) x, @sigmoid, @exp, @scaledSigmoid, @sigmoid};

param_names{10} = {'\epsilon', '\rho', 'goBias','\beta_{\Omega}', '\thres_{\Omega}', 'reward info'}; % DynamicOmega2Model
transform{10} = {@sigmoid, @exp, @(x) x, @sigmoid, @exp, @scaledSigmoid, @(x) x};

param_names{11} = {'\epsilon', '\rho', 'goBias'};
transform{11} = {@sigmoid, @exp, @(x) x};

param_names{12} = {'\epsilon', '\rho', 'goBias', '\alpha_{up}', '\alpha_{down}'};
transform{12} = {@sigmoid, @exp, @(x) x, @sigmoid, @sigmoid};

param_names{13} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{13} = {@sigmoid, @exp, @(x) x, @sigmoid, @exp, @scaledSigmoid, @sigmoid};

% ----------------------------------------------------------------------- %
%% 00b) Select model:

% Get list of .m files in dir.models
modelFiles = dir(fullfile(dirs.models, '*.m'));

% Set nMod to the number of .m files
nMod = length(modelFiles);

selMod = 3;

fprintf("Selected model is no. %02d\n", selMod)

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 01) LOAD AND SAVE:
%% 01a) Alternative A: Load model fitted with LAP:

parType = 'lap';

fname_mod = cell(nMod, 1);

for iMod = 1:nMod
    fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d.mat', iMod));
end

fprintf("Load model %d fit with LAP\n", selMod);
fname       = load(fname_mod{selMod});
cbm         = fname.cbm;
subParam    = cbm.output.parameters;
nSub        = size(subParam, 1);
nParam      = size(subParam, 2);

% ----------------------------------------------------------------------- %
%% 01b) Alternative B: Load model fitted with HBI:

parType = 'hbi';

% load
fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d.mat', selMod));
fprintf("Load model %d fit with HBI\n", selMod);
load(fname_hbi)

% extract
fprintf("Extract model %d fit with HBI\n", selMod);
groupParam      = cbm.output.group_mean{:};
subParam        = cbm.output.parameters{:};
nSub            = size(subParam, 1);
nParam          = size(subParam, 2);



% ----------------------------------------------------------------------- %
%% 01c) Transform group and subject parameters appropriately:
% Define sigmoid function
sigmoid = @(x) 1 ./ (1 + exp(-x));

fprintf('Transform parameters of model %d\n', selMod);

for iParam = 1:nParam
    % Get the transformation function for the current parameter
    transformFunc = transform{selMod}{iParam};
    groupParam(iParam) = transformFunc(groupParam(iParam));
    subParam(:, iParam) = transformFunc(subParam(:, iParam));
end

% ----------------------------------------------------------------------- %
%% 01d) Save (either LAP or HBI):

fprintf('Save transformed group and subject level parameters under %s, model %02d\n', ...
    parType, selMod);

fileName = fullfile(eval(sprintf('dirs.%s', parType)), sprintf('CBM_%s_M%02d_groupParameters.csv', ...
    parType, selMod));

csvwrite(fileName, groupParam);

fileName = fullfile(eval(sprintf('dirs.%s', parType)), sprintf('CBM_%s_M%02d_subjectParameters.csv', ...
    parType, selMod));

csvwrite(fileName, subParam);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 02) GROUP-LEVEL PARAMETERS.
%% 02a) Mean per proup level parameter:

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
%% 03) SUBJECT-LEVEL PARAMETERS.
%% 03a) Mean, SD, range per subject-level parameter:

fprintf('\nSubject-level parameters for model %d\n', selMod);
for iParam = 1:size(subParam, 2)
    x = subParam(:, iParam);
    fprintf('Parameter %d: M = %.03f, SD = %.03f, range %.03f - %.03f \n', ...
        iParam, mean(x), std(x), min(x), max(x));
end

% ----------------------------------------------------------------------- %
%% 03b) Find extreme values:

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
%% 03c) Sign of parameter:

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
%% 04) PLOTS.
%% 04a) Plot parameters with cbm function:
% Parameter names:
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

% Model names:
model_names = {'M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12', 'M13'};
% Create output name:
modVec = 1:nMod;
fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s.mat', num2str(modVec, '_%02d')));
d = load(fname_hbi);
c = d.cbm;
freq = c.output.model_frequency;
[~, k] = max(freq);

cbm_hbi_plot(fname_hbi, model_names, param_names{k}, transform{k});

% ----------------------------------------------------------------------- %
%% 04b) Density plots of single parameters:

f_hbi   = load(fname_hbi);
cbm     = f_hbi.cbm;
selMod    = 5; % specify model to look at
% Click through densityplots of each parameter:
for iParam = 1:size(cbm.output.parameters{selMod}, 2) % iParam = 2
    x = cbm.output.parameters{selMod}(:, iParam);
    ksdensity(x); hold on
    plot(x, 0.1 + randn(length(x), 1)/50, 'b.');
    title(sprintf('Model %02d: Parameter %d', selMod, iParam));
    w = waitforbuttonpress;
    close gcf
end

% ----------------------------------------------------------------------- %
%% 04c) Correlations/ bivariate distributions:

% Select parameters:
iParam1 = 1; % specify first parameter
iParam2 = 4; % specify second parameter
% Retrieve parameter values:
x       = subParam(:, iParam1);
y       = subParam(:, iParam2);
p       = polyfit(x', y', 1);
yhat    = polyval(p, x);

% Plot:
plot(x, y, 'b.'); hold on;
plot(x, yhat, 'r-');
axis([min(x) max(x) min(y) max(y)]);

% Print correlation to console:
corVal  = corr(x, y);
fprintf('Model %02d: Parameters %d and %d, cor = %0.2f\n', ...
    selMod, iParam1, iParam2, corVal);
title(sprintf('Model %02d: Parameters %d and %d, cor = %0.2f', ...
    selMod, iParam1, iParam2, corVal));

% ----------------------------------------------------------------------- %
%% 04d) Bar plots with dots:

% Settings:
paramNames  = {{'$\rho$', '$\epsilon$'},...
    {'$\rho$', '$\epsilon$', '{\it Go}'},...
    {'$\rho$', '$\epsilon$', '{\it Go}', '$\pi$'},...
    {'$\rho$', '$\epsilon$', '{\it Go}', '$\kappa$'},...
    {'$\rho$', '$\epsilon$', '{\it Go}', '$\pi$', '$\kappa$'}};
ScatterMatrix = subParam';
colMat = [0 0 1; 1 0 0; 0 01 0; 0 1 1; 1 0 1;1 0 0; 1 0 0]; % color
posMat = 1:1:(0.5 + nParam);

% Make figure:
addpath(fullfile(dirs.root, '/Analyses/Behavior_Scripts/Matlab_Plots/')); % add function barScatter
figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold on
barScatter(ScatterMatrix, [], [], true, true, colMat, posMat);

% Add plot features:
set(gca, 'xlim', [0 (nParam+1)], 'ylim', [-5 5],...
    'xtick', 1:nParam, 'xticklabel', paramNames{selMod},...
    'FontSize', 32, 'FontName', 'Arial', 'Linewidth', 4, 'TickLabelInterpreter', 'latex');
xlabel('Parameter', 'FontSize',32, 'FontName', 'Arial');
ylabel('Parameter estimates', 'FontSize',32, 'FontName', 'Arial');
box off; hold off

% Print to console:
fprintf('Model %02d:\n', selMod)
for ii = 1:size(ScatterMatrix, 1)
    paramVec = ScatterMatrix(ii, :);
    fprintf('Parameter %02d: M = %.02f, SD = %.02f, range %.02f - %.02f\n', ...
        ii, nanmean(paramVec), nanstd(paramVec), min(paramVec), max(paramVec));
end

% Save:
saveas(gcf, fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/Parameters_%s_mod%02d.png', ...
    parType, selMod)));
pause(2)
close gcf

% END OF FILE.