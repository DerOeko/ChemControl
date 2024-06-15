% ChemControl_cbm_eval_priors

% This is an interactive script---execute it step-by-step.m
% - First: In which ranges can our recovery process not recover parameters.
% 	- For this, we need matching priors for the recovery and sampling process. We just want to know, "If we know the exact same priors, can we recover the parameters?".
% - The goal of this: Find the least constraining priors, for which we can still recover parameters.
% - Second: In which ranges can our model recovery process not recover models?
%	- For this, we need the priors from the previous step, and try to retrieve the correct model. The further away the resulting confusion matrix is from the identity matrix, the more we have to constrain the priors.
% - Third: Take the least constraining priors for the parameter and model recovery that still gave good results. Fit the real data using those priors. Plot the correlation between the simulated learning curves of the winning model and the learning curves of the average participant.
% - Fourth: If they match, and the winning model is omega, stop and be happy.
% - If they match, and the winning model isn't omega, find a new model.
% - If they don't match, return to improving priors.
%
% EXPLANATION OF SETTINGS:
% dirs.root               = string, root directory of project.
%
% OUTPUTS:
% no outputs, just plots.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel
clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Set directories:
dirs = [];
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';

fprintf('Initialize directories\n');

dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

dirs.lap     = fullfile(dirs.results, 'LAP_Results');
if ~exist(dirs.lap, 'dir'); mkdir(dirs.lap); end

dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
if ~exist(dirs.hbi, 'dir'); mkdir(dirs.hbi); end

dirs.models         = fullfile(dirs.root, 'models');

% Create new directories if necessary:
dirs.paramRecovery = fullfile(dirs.results, 'parameterRecovery');
if ~exist(dirs.paramRecovery, 'dir'); mkdir(dirs.paramRecovery); end


%% 00b) Settings:
% Get list of .m files in dir.models
modelFiles = dir(fullfile(dirs.models, '*.m'));

% Set nMod to the number of .m files
nMod = length(modelFiles);

nParams = [2 3 4 4 4 7 7]; % Number of parameters per model
selMod = 7; % Which model?
nParam = nParams(selMod); % Number of params for this model

numSamples = 100; % How many parameter samples for each parameter?
numSampleParam = 100; % From each parameter samples, how many to take
nSub = 1; % How many simulations for each parameter combination?

modelSimHandle = str2func(sprintf('ChemControl_mod%d_modSim', selMod));
modelHandle = str2func(sprintf('ChemControl_mod%d', selMod));
fprintf('Selected model %d and number of parameters %d\n', selMod, nParam)

%% 01a) Prior finding with Grid Search
% Define ranges for grid search
mean_range = -2:0.5:2; % Range of means for grid search
variance_range = 1:1:10; % Range of variances for grid search

eps_mean    = 0; eps_v      = 5; % Learning rate (sigmoid)
rho_mean    = 2; rho_v      = 3; % Discount factor (exp)
gB_mean     = 0; gB_v       = 10; % Go bias (identity)
slope_mean  = 2; slope_v    = 3; % Slope (exp)

% Initialize result storage
result_rows = length(mean_range) * length(variance_range);
recoverability_scores = zeros(result_rows, 3);

% Initialize row counter
row_counter = 1;
for mu = 1:length(mean_range)
    for var = 1:length(variance_range)
        % Set current priors
        current_mean = mean_range(mu);
        current_var = variance_range(var);
        
        % Update priors for the current combination
        priors{6}   = struct('mean', [eps_mean, rho_mean, gB_mean, current_mean, current_mean, current_mean, slope_mean], ...
                            'variance', [eps_v, rho_v, gB_v, current_var, current_var, current_var, slope_v]); % prior_model_dynamicOmega1
        eps     = normrnd(eps_mean,   sqrt(eps_v), [numSamples, 1]);
        rhos    = normrnd(rho_mean,   sqrt(rho_v), [numSamples, 1]);
        gBs     = normrnd(gB_mean,    sqrt(gB_v), [numSamples, 1]);
        o0s     = normrnd(current_mean,    sqrt(current_var), [numSamples, 1]);
        alphas  = normrnd(current_mean, sqrt(current_var), [numSamples, 1]);
        kaps    = normrnd(current_mean,   sqrt(current_var), [numSamples, 1]);
        slopes  = normrnd(slope_mean, sqrt(slope_v), [numSamples, 1]);
        
        % Create parameter combinations
        paramCombinations = [randsample(eps, numSampleParam, true), ...
                             randsample(rhos, numSampleParam, true), ...
                             randsample(gBs, numSampleParam, true), ...
                             randsample(o0s, numSampleParam, true), ...
                             randsample(alphas, numSampleParam, true), ...
                             randsample(kaps, numSampleParam, true), ...
                             randsample(slopes, numSampleParam, true)];

        numCombinations = uint64(size(paramCombinations, 1));
        % Run parameter recovery
        retrievedParams = zeros(numCombinations, nParam);

        parfor iComb = 1:numCombinations
            parameters = paramCombinations(iComb, 1:nParam);
            data = cell(nSub, 1);
        
            for iSub = 1:nSub
                % Extract subject data:
                fprintf("Start subject %03d\n", iSub)
                subj = sim_subj;
                % Simulate:
                out = modelSimHandle(parameters, subj); 
                data{iSub} = struct("stimuli", out.stimuli, "actions", out.actions, "outcomes", out.outcomes);
            end
        
            temp_folder = tempname;
            mkdir(temp_folder);
            fname = fullfile(temp_folder, 'lap.mat');
        
            cbm_lap(data, modelHandle, priors{selMod}, fname);
            d = load(fname);
            params = d.cbm.output.parameters;
        
            % Store temporary results into the final matrix
            retrievedParams(iComb, :) = params(:);
        end

        % Calculate recoverability score (correlation)
        corr_values = zeros(1, nParam);
        for iParam = 1:nParam
            corr_values(iParam) = corr(paramCombinations(:, iParam), retrievedParams(:, iParam));
        end
        avg_correlation = mean(corr_values);

        % Store the current priors and recoverability score
        recoverability_scores(row_counter, :) = [current_mean, current_var, avg_correlation];
        row_counter = row_counter + 1;
    end
end
% Save results to a file
save(fullfile(dirs.results, 'recoverability_scores.mat'), 'recoverability_scores');

%% 01b) Summary statistics
T = load(fullfile(dirs.results, 'recoverability_scores.mat'));
T = T.recoverability_scores;
low_recovery = T(T(:, 3) <= 0.5, :);

summary_stats_mean = array2table(mean(low_recovery), 'VariableNames', {'mean', 'var', 'score'});
disp('Mean values for parameters with recoverability score below 0.5:');
disp(summary_stats_mean);

summary_stats_median = array2table(median(low_recovery), 'VariableNames', {'mean', 'var', 'score'});
disp('Median values for parameters with recoverability score below 0.5:');
disp(summary_stats_median);


summary_stats_mean = array2table(mean(T), 'VariableNames', {'mean', 'var', 'score'});
disp('Mean values for parameters with recoverability score above 0.5:');
disp(summary_stats_mean);

summary_stats_median = array2table(median(T), 'VariableNames', {'mean', 'var', 'score'});
disp('Median values for parameters with recoverability score above 0.5:');
disp(summary_stats_median);

% ----------------------------------------------------------------------- %
%% 01c) Select good parameters and simulate recovery
% Settings:
fprintf('>>> Set parameter recovery correlation settings\n')
numSamples = 1000; % How many parameter samples for each parameter?
numSampleParam = 10000; % From each parameter samples, how many to take

% Define means and variances for each parameter
eps_mean    = 0; eps_v      = 10; % Learning rate (sigmoid)
rho_mean    = 2; rho_v      = 3; % Discount factor (exp)
gB_mean     = 0; gB_v       = 10; % Go bias (identity)
pi_mean     = 0; pi_v       = 10; % Pavlovian bias (identity)
o_mean      = 0; o_v        = 5; % Omega (sigmoid)
o0_mean     = 5; o0_v       = 100; % Initial omega (sigmoid)
alpha_mean  = 0; alpha_v    = 5; % Alpha (sigmoid)
kap_mean    = 0; kap_v      = 5; % Kappa (sigmoid)
slope_mean  = 2; slope_v    = 3; % Slope (exp)
aO_mean     = 0; aO_v       = 100; % Alpha Omega (sigmoid)
bO_mean     = 2; bO_v       = 100; % Beta Omega (exp)
tO_mean     = 0; tO_v       = 100; % Theta Omega (scaled sigmoid)

% Generate samples for each parameter
eps     = normrnd(eps_mean,   sqrt(eps_v), [numSamples, 1]);
rhos    = normrnd(rho_mean,   sqrt(rho_v), [numSamples, 1]);
gBs     = normrnd(gB_mean,    sqrt(gB_v), [numSamples, 1]);
pis     = normrnd(pi_mean,    sqrt(pi_v), [numSamples, 1]);
os      = normrnd(o_mean,     sqrt(o_v), [numSamples, 1]);
o0s     = normrnd(o0_mean,    sqrt(o0_v), [numSamples, 1]);
alphas  = normrnd(alpha_mean, sqrt(alpha_v), [numSamples, 1]);
kaps    = normrnd(kap_mean,   sqrt(kap_v), [numSamples, 1]);
slopes  = normrnd(slope_mean, sqrt(slope_v), [numSamples, 1]);
aOs     = normrnd(aO_mean,    sqrt(aO_v), [numSamples, 1]);
bOs     = normrnd(bO_mean,    sqrt(bO_v), [numSamples, 1]);
tOs     = normrnd(tO_mean,    sqrt(tO_v), [numSamples, 1]);

% Create parameter combinations
paramCombinations = [randsample(eps, numSampleParam, true), ...
                     randsample(rhos, numSampleParam, true), ...
                     randsample(gBs, numSampleParam, true), ...
                     randsample(pis, numSampleParam, true), ...
                     randsample(os, numSampleParam, true), ...
                     randsample(o0s, numSampleParam, true), ...
                     randsample(alphas, numSampleParam, true), ...
                     randsample(kaps, numSampleParam, true), ...
                     randsample(slopes, numSampleParam, true), ...
                     randsample(aOs, numSampleParam, true), ...
                     randsample(bOs, numSampleParam, true), ...
                     randsample(tOs, numSampleParam, true)];

% paramCombinations = [randsample(eps, numSampleParam, true), ...
%                      randsample(rhos, numSampleParam, true), ...
%                      randsample(gBs, numSampleParam, true), ...
%                      randsample(o0s, numSampleParam, true), ...
%                      randsample(aOs, numSampleParam, true), ...
%                      randsample(bOs, numSampleParam, true), ...
%                      randsample(tOs, numSampleParam, true)];

numCombinations     = uint64(size(paramCombinations, 1));

% Transform functions
transform{1} = {@sigmoid, @exp};
transform{2} = {@sigmoid, @exp, @(x) x};
transform{3} = {@sigmoid, @exp, @(x) x, @(x) x};
transform{4} = {@sigmoid, @exp, @(x) x, @(x) x};
transform{5} = {@sigmoid, @exp, @(x) x, @sigmoid};
transform{6} = {@sigmoid, @exp, @(x) x, @sigmoid, @sigmoid, @sigmoid, @exp};
transform{7} = {@sigmoid, @exp, @(x) x, @sigmoid, @sigmoid, @exp, @scaledSigmoid};

% Priors:
fprintf('>>> Initialize unconstraining parameter priors\n')
priors{1}   = struct('mean', [eps_mean, rho_mean], ...
                    'variance', [eps_v, rho_v]); % prior_model

priors{2}   = struct('mean', [eps_mean, rho_mean, gB_mean], ...
                    'variance', [eps_v, rho_v, gB_v]); % prior_model_goBias

priors{3}   = struct('mean', [eps_mean, rho_mean, gB_mean, pi_mean], ...
                    'variance', [eps_v, rho_v, gB_v, pi_v]); % prior_model_fixedPavlov

priors{4}   = struct('mean', [eps_mean, rho_mean, gB_mean, pi_mean], ...
                    'variance', [eps_v, rho_v, gB_v, pi_v]); % prior_model_dynamicPavlov

priors{5}   = struct('mean', [eps_mean, rho_mean, gB_mean, o_mean], ...
                    'variance', [eps_v, rho_v, gB_v, o_v]); % prior_model_fixedOmega

priors{6}   = struct('mean', [eps_mean, rho_mean, gB_mean, o0_mean, alpha_mean, kap_mean, slope_mean], ...
                    'variance', [eps_v, rho_v, gB_v, o0_v, alpha_v, kap_v, slope_v]); % prior_model_dynamicOmega1

priors{7}   = struct('mean', [eps_mean, rho_mean, gB_mean, o0_mean, aO_mean, bO_mean, tO_mean], ...
                    'variance', [eps_v, rho_v, gB_v, o0_v, aO_v, bO_v, tO_v]); % prior_model_dynamicOmega2

% Define parameter indices for each model
paramIndices = { [1, 2], ...
                 [1, 2, 3], ...
                 [1, 2, 3, 4], ...
                 [1, 2, 3, 4], ...
                 [1, 2, 3, 5], ...
                 [1, 2, 3, 6, 7, 8, 9], ...
                 [1, 2, 3, 6, 10, 11, 12]};

% Select the correct parameter indices for the current model
currentParamIndices = paramIndices{selMod};

% Run parameter recovery
retrievedParams = zeros(numCombinations, nParam);

parfor iComb = 1:numCombinations
    % parameters = paramCombinations(iComb, :);
    parameters = zeros(nParam, 1);
    for i = 1:nParam
        iParam = currentParamIndices(i);
        parameters(i) = paramCombinations(iComb, iParam);
    end
    data = cell(nSub, 1);

    for iSub = 1:nSub
        % Extract subject data:
        fprintf("Start subject %03d\n", iSub)
        subj = sim_subj;
        % Simulate:
        out = modelSimHandle(parameters, subj); 
        data{iSub} = struct("stimuli", out.stimuli, "actions", out.actions, "outcomes", out.outcomes);
    end

    temp_folder = tempname;
    mkdir(temp_folder);
    fname = fullfile(temp_folder, 'lap.mat');

    cbm_lap(data, modelHandle, priors{selMod}, fname);
    d = load(fname);
    params = d.cbm.output.parameters;

    % Store temporary results into the final matrix
    retrievedParams(iComb, :) = params(:);
end

%% 01d) Plotting of results with correlation lines
figure;
sgtitle('Parameter Recovery Analysis');


% Settings:
paramNames = {{'\epsilon', '\rho'},...
    {'\epsilon', '\rho', 'goBias'},...
    {'\epsilon', '\rho', 'goBias', '\pi'},...
    {'\epsilon', '\rho', 'goBias', '\pi'},...
    {'\epsilon', '\rho', 'goBias', '\omega'},...
    {'\epsilon', '\rho', 'goBias', '\omega_{init}', '\alpha', '\kappa', 'slope'},...
    {'\epsilon', '\rho', 'goBias', '\omega_{init}', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}'}};

for i = 1:nParam
    iParam = currentParamIndices(i);
    disp(iParam)
    subplot(1, nParam, i);
    scatter(paramCombinations(:, iParam), retrievedParams(:, i), "filled");
    hold on;
    % Fit a linear model to the data
    p = polyfit(paramCombinations(:, iParam), retrievedParams(:, i), 1);
    yfit = polyval(p, paramCombinations(:, iParam));
    % Plot the linear fit
    plot(paramCombinations(:, iParam), yfit, 'LineWidth', 2);
    xlabel(sprintf('True %s', paramNames{selMod}{i}));
    ylabel(sprintf('Retrieved %s', paramNames{selMod}{i}));
    title(sprintf('True vs. Retrieved %s, \ncorr: %.02f', paramNames{selMod}{i}, corr(paramCombinations(:, iParam), retrievedParams(:, i))));
end
% for i = 1:nParam
%     subplot(1, nParam, i);
%     scatter(paramCombinations(:, i), retrievedParams(:, i), "filled");
%     hold on;
%     % Fit a linear model to the data
%     p = polyfit(paramCombinations(:, i), retrievedParams(:, i), 1);
%     yfit = polyval(p, paramCombinations(:, i));
%     % Plot the linear fit
%     plot(paramCombinations(:, i), yfit, 'LineWidth', 2);
%     xlabel(sprintf('True %s', paramNames{selMod}{i}));
%     ylabel(sprintf('Retrieved %s', paramNames{selMod}{i}));
%     title(sprintf('True vs. Retrieved %s, \ncorr: %.02f', paramNames{selMod}{i}, corr(paramCombinations(:, i), retrievedParams(:, i))));
% end

% ----------------------------------------------------------------------- %
%% 01e) Correlations/ bivariate distributions:

% Select parameters:
iParam1 = 1; % specify first parameter
iParam2 = 2; % specify second parameter

for iParam1 = 1:nParam
    for iParam2 = 1:nParam
        if iParam1 ~= iParam2
            % Retrieve parameter values:
            x       = retrievedParams(:, iParam1);
            y       = retrievedParams(:, iParam2);
            p       = polyfit(x', y', 1);
            yhat    = polyval(p, x);
            
            % Plot:
            plot(x, y, 'b.'); hold on;
            plot(x, yhat, 'r-');
            axis([min(x) max(x) min(y) max(y)]);
            xlabel(paramNames{selMod}{iParam1});
            ylabel(paramNames{selMod}{iParam2});
            
            % Print correlation to console:
            corVal  = corr(x, y);
            fprintf('Model %02d: Parameters %d and %d, cor = %0.2f\n', ...
            selMod, iParam1, iParam2, corVal);
        title(sprintf('Model %02d: Parameters %d and %d, cor = %0.2f', ...
            selMod, iParam1, iParam2, corVal));
            w = waitforbuttonpress;
            close gcf
        end
    end
end

% ----------------------------------------------------------------------- %
%% 02a) Confusion matrix: 

confusionMatrix = zeros(nMod, nMod);
runsPerM = 10;
Mods = repelem([1 2 3 4 5 6 7], runsPerM);
nRuns = nMod * runsPerM;

% Usa a temporary cell array to store intermediate confusion matrices
% (needed for parallel computing)
tempConfusionMatrices = cell(nRuns, 1);

for iMod = 1:nMod
    nParam = nParams(iMod);
    currentParamIndices = paramIndices(iMod);
    parfor iRun = 1:runsPerM
        parameters = zeros(nParam, 1);
        for i = 1:nParam
            iParam = currentParamIndices(i);
            parameters(i) = paramCombinations(iComb, iParam);
        end
        data = cell(nSub, 1);
        for iSub = 1:nSub
            % Extract subject data:
            fprintf("Start subject %03d\n", iSub)
            subj = sim_subj;
            % Simulate:
            out = modelSimHandle(parameters, subj); 
            data{iSub} = struct("stimuli", out.stimuli, "actions", out.actions, "outcomes", out.outcomes);
        end

        temp_folder = tempname;
        mkdir(temp_folder);
        fcbm_maps = cell(nMod, 1);
        models = cell(nMod, 1);
        for iMod = 1:nMod
            fcbm_maps{iMod} = fullfile(temp_folder, sprintf("lap_mod%02d.mat", iMod));
            models{iMod} = str2func(sprintf("@ChemControl_mod%d", iMod));
            cbm_lap(data, models{iMod}, priors{iMod}, fcbm_maps{iMod});
        end

        fname_hbi = fullfile(temp_folder, 'hbi_model.mat');
        cbm_hbi(data, models, fcbm_maps, fname_hbi)

        d = load(fname_hbi);
        cbm = d.cbm;
        freq = cbm.output.model_frequency;
    
        [~, choice] = max(freq);
    
        % Store the result in the temporary confusion matrix
        tempConfusionMatrices{iRun} = zeros(nMod, nMod);
        tempConfusionMatrices{iRun}(iMod, choice) = 1;
    end
end

for iRun = 1:nRuns
    confusionMatrix = confusionMatrix + tempConfusionMatrices{iRun};
end

% Display the confusion matrix
disp('Confusion Matrix:')
disp(confusionMatrix)

% Display confusion matrix as a heatmap
figure;

% Create the heatmap
h = heatmap(confusionMatrix, 'Colormap', parula, 'ColorbarVisible', 'on', 'Title', 'Confusion Matrix');
h.XDisplayLabels = cfg.model_names; % Set x-axis labels
h.YDisplayLabels = cfg.model_names; % Set y-axis labels

% Add labels to the axes
ax = gca;

% Customize the heatmap axes to place x-axis labels at the top
h.NodeChildren(3).XAxisLocation = 'top';
xlabel('Predicted Model');
ylabel('True Model');

        