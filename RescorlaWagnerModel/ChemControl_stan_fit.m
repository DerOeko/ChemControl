% This is an interactive script---execute it step-by-step.
% It fits a series of computational reinforcement learning models using the
% CBM toolbox and evaluates them.
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
%
% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Directories:
dirs = [];
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';

fprintf('Initialize directories\n');

dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

dirs.lap     = fullfile(dirs.results, 'LAP_Results');
if ~exist(dirs.lap, 'dir'); mkdir(dirs.lap); end

dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
if ~exist(dirs.hbi, 'dir'); mkdir(dirs.hbi); end

dirs.stan = fullfile(dirs.results, 'Stan_Results');
if ~exist(dirs.stan, 'dir'); mkdir(dirs.stan); end

dirs.models         = fullfile(dirs.root, 'models');

dirs.stan_models = fullfile(dirs.root, 'stan_models');
if ~exist(dirs.stan_models, 'dir'); mkdir(dirs.stan_models); end

% ----------------------------------------------------------------------- %
%% 00b) Settings:
nMod = length(dir(fullfile(dirs.models, "*.m")));

fprintf('Fit %d models\n', nMod);

% ----------------------------------------------------------------------- %
%% 00c) Paths:
fprintf('Add paths\n');

% Add paths: 
addpath('/home/control/samnel/Documents/MATLAB/cbm-master/codes'); % CBM toolbox
addpath(fullfile(dirs.root, 'behavioral_study', 'scripts', 'matlab_scripts', 'RescorlaWagnerModel', 'models')); % models

% ----------------------------------------------------------------------- %
%% 00d) Simulation settings:

nParams = [2 3 4 4 4 7 7]; % Number of parameters per model
selMod = 7; % Which model?
nParam = nParams(selMod); % Number of params for this model

numSamples = 100; % How many parameter samples for each parameter?
numSampleParam = 100; % From each parameter samples, how many to take
nSub = 1; % How many simulations for each parameter combination?

modelSimHandle = str2func(sprintf('ChemControl_mod%d_modSim', selMod));
modelHandle = str2func(sprintf('ChemControl_mod%d', selMod));
fprintf('Selected model %d and number of parameters %d\n', selMod, nParam)

% Define means and variances for each parameter
eps_mean    = 0; eps_v      = 2; % Learning rate (sigmoid)
rho_mean    = 1; rho_v      = 1; % Discount factor (exp)
gB_mean     = 0; gB_v       = 2; % Go bias (identity)
oi_mean     = 0; oi_v       = 2; % Omega init value (sigmoid)
aO_mean     = 0; aO_v       = 0.5; % Alpha Omega (sigmoid)
bO_mean     = 1; bO_v       = 1; % Beta Omega (exp)
tO_mean     = 0; tO_v       = 2; % Theta Omega (scaled sigmoid)

% Generate samples for each parameter
eps     = sigmoid(normrnd(eps_mean,   sqrt(eps_v), [numSamples, 1]));
rhos    = exp(normrnd(rho_mean,   sqrt(rho_v), [numSamples, 1]));
gBs     = normrnd(gB_mean,    sqrt(gB_v), [numSamples, 1]);
ois     = sigmoid(normrnd(oi_mean, sqrt(oi_v), [numSamples, 1]));
aOs     = sigmoid(normrnd(aO_mean,    sqrt(aO_v), [numSamples, 1]));
bOs     = exp(normrnd(bO_mean,    sqrt(bO_v), [numSamples, 1]));
tOs     = scaledSigmoid(normrnd(tO_mean,    sqrt(tO_v), [numSamples, 1]));

% Create parameter combinations
paramCombinations = [randsample(eps, numSampleParam, true), ...
                     randsample(rhos, numSampleParam, true), ...
                     randsample(gBs, numSampleParam, true), ...
                     randsample(ois, numSampleParam, true), ...
                     randsample(aOs, numSampleParam, true), ...
                     randsample(bOs, numSampleParam, true), ...
                     randsample(tOs, numSampleParam, true)];

numCombinations     = uint64(size(paramCombinations, 1));

nRuns = 100;

%% 00e) Fitting settings:
numChains = maxNumCompThreads;
numIter = 2000;

%% 00f) Data and parameters for fitting
nBlocks = 8;
nTrials = 40;
nSubjects = 1;
retrievedParams = zeros(nRuns, nParam);
trueParams = zeros(nRuns, nParam);
for iRun = 1:nRuns
    fprintf('Fit %03d started\n', iRun);
    subj = sim_subj();
    parameters = paramCombinations(iRun, :);

    out = modelSimHandle(parameters, subj);
    states = reshape(out.stimuli, [nSubjects, nBlocks, nTrials]);
    actions = reshape(out.actions, [nSubjects, nBlocks, nTrials]);
    actions(actions == 2) = 0;
    feedbacks = reshape(out.outcomes, [nSubjects, nBlocks, nTrials]);
    d = struct('nTrials', nTrials, 'nBlocks', nBlocks, 'nSubjects', nSubjects, 'state', states, 'action', actions, 'feedback', feedbacks);
    fit = stan('file', fullfile(dirs.stan_models, sprintf('ChemControl_stan_mod%02d.stan', selMod)), 'working_dir', '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/Log/Behavior/Modelling_CBM/Stan_Results/', 'verbose', false, 'data', d, 'chains', numChains, 'iter', numIter);
    fit.block();
    fprintf('\n---------------------------------------------------\n');
    fprintf('Fit complete\n');

    % 01a) Extract parameters:
    paramNames = {'ep', 'rho', 'goBias', 'omegaInit', 'alphaO', 'betaO', 'thresO'};
    samples = fit.extract('par', paramNames);
    epsilon_samples = reshape(samples.ep, [numIter*nBlocks/2 * numChains, 1]);
    rho_samples = reshape(samples.rho, [numIter*nBlocks/2 * numChains, 1]);
    goBias_samples = reshape(samples.goBias, [numIter*nBlocks/2 * numChains, 1]);
    omegaInit_samples = reshape(samples.omegaInit, [numIter*nBlocks/2 * numChains, 1]);
    alphaO_samples = reshape(samples.alphaO, [numIter*nBlocks/2 * numChains, 1]);
    betaO_samples = reshape(samples.betaO, [numIter*nBlocks/2 * numChains, 1]);
    thresO_samples = reshape(samples.thresO, [numIter*nBlocks/2 * numChains, 1]);
    
    epsilon_median = median(epsilon_samples);
    rho_median = median(rho_samples);
    goBias_median = median(goBias_samples);
    omegaInit_median = median(omegaInit_samples);
    alphaO_median = median(alphaO_samples);
    betaO_median = median(betaO_samples);
    thresO_median = median(thresO_samples);
    
    retrievedParams(iRun, :) = [epsilon_median, rho_median, goBias_median, omegaInit_median, alphaO_median, betaO_median, thresO_median];
    trueParams(iRun, :) = parameters;
end

% 01b) Plotting of results with correlation lines
figure;
sgtitle('Parameter Recovery Analysis');


% Settings:
paramNames = {{'\epsilon', '\rho', 'goBias', '\omega_{init}', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}'}};

for i = 1:nParam

    subplot(1, nParam, i);
    scatter(trueParams(:, i), retrievedParams(:, i), "filled");
    hold on;
    % Fit a linear model to the data
    p = polyfit(trueParams(:, i), retrievedParams(:, i), 1);
    yfit = polyval(p, trueParams(:, i));
    % Plot the linear fit
    plot(trueParams(:, i), yfit, 'LineWidth', 2);

    % Plot the optimal line (y = x)
    xLimits = xlim;
    yLimits = ylim;
    plot(xLimits, xLimits, 'k--', 'LineWidth', 2); % Diagonal optimal line

    xlabel(sprintf('True %s', paramNames{1}{i}));
    ylabel(sprintf('Retrieved %s', paramNames{1}{i}));
    title(sprintf('True vs. Retrieved %s, \ncorr: %.02f', paramNames{1}{i}, corr(trueParams(:, i), retrievedParams(:, i))));
end

% True parameter values
% epsilon_true = 0.3;
% rho_true = 5;
% goBias_true = 0.5;
% omegaInit_true = 0.3;
% alphaO_true = 0.1;
% betaO_true = 5;
% thresO_true = 0.5;
% 
% Fancy print statements
% fprintf('\n---------------------------------------------------\n');
% fprintf('Comparison of estimated and true parameter values:\n');
% fprintf('---------------------------------------------------\n');
% fprintf('Parameter\t\tEstimated\t\tTrue\n');
% fprintf('---------------------------------------------------\n');
% fprintf('Epsilon\t\t%.4f\t\t%.4f\n', epsilon_median, epsilon_true);
% fprintf('Rho\t\t\t%.4f\t\t%.4f\n', rho_median, rho_true);
% fprintf('GoBias\t\t\t%.4f\t\t%.4f\n', goBias_median, goBias_true);
% fprintf('omegaInit\t\t%.4f\t\t%.4f\n', omegaInit_median, omegaInit_true);
% fprintf('alphaO\t\t\t%.4f\t\t%.4f\n', alphaO_median, alphaO_true);
% fprintf('betaO\t\t\t%.4f\t\t%.4f\n', betaO_median, betaO_true);
% fprintf('thresO\t\t\t%.4f\t\t%.4f\n', thresO_median, thresO_true);
% fprintf('---------------------------------------------------\n');
% 
% figure;
% hold on;
% histogram(epsilon_samples, 150, 'DisplayName', 'Epsilon Samples');
% histogram(rho_samples, 150, 'DisplayName', 'Rho Samples');
% histogram(goBias_samples, 150, 'DisplayName', 'GoBias Samples');
% histogram(omegaInit_samples, 150, 'DisplayName', 'OmegaInit Samples');
% histogram(alphaO_samples, 150, 'DisplayName', 'AlphaO Samples');
% histogram(betaO_samples, 150, 'DisplayName', 'BetaO Samples');
% 
% hold off
% xlabel('Parameter estimate')
% ylabel('Count')
% title('Parameter retrieval results using MCMC sampling')
% legend('Location', 'best')
% 
% xline(epsilon_true, 'LineWidth', 1, 'LineStyle', ':', 'Color', 'r', 'DisplayName', 'True Epsilon');
% xline(rho_true, 'LineWidth', 1, 'LineStyle', ':', 'Color', 'b', 'DisplayName', 'True Rho');
% xline(goBias_true, 'LineWidth', 1, 'LineStyle', ':', 'Color', 'g', 'DisplayName', 'True GoBias');
% xline(omegaInit_true, 'LineWidth', 1, 'LineStyle', ':', 'Color', 'm', 'DisplayName', 'True OmegaInit');
% xline(alphaO_true, 'LineWidth', 1, 'LineStyle', ':', 'Color', 'y', 'DisplayName', 'True AlphaO');
% xline(betaO_true, 'LineWidth', 1, 'LineStyle', ':', 'Color', 'c', 'DisplayName', 'True BetaO');
% 
% legend('Location', 'best');