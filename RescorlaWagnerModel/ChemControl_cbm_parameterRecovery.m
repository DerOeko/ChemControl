% ChemControl_cbm_parameterRecovery.m

% This is an interactive script---execute it step-by-step.
% Its performs parameter recovery for a selected model.
% The following steps have to be executed beforehand:
% - Fit model using ChemControl_cbm_fit to obtain parameters fitted to
% actual data.
% - Simulate new data using ChemControl_cbm_sim given parameters fitted.
% Mind to set the root directory dirs.root.

%
% EXPLANATION OF SETTINGS:
% dirs.root     = string, root directory of project.
% selMod          = numeric integer, ID of selected model for which to
% perform model recovery.
% simType       = scalar string, type of simulations, should be 'modSim'
% for probabilistic simulations of actual responses. 'realSim' if model sim
% should be based on real data.
% parType       = scalar string, type of parameters used for simulations,
% either LaPlace approximation ('LAP') or hierarchical Bayesian inference
% (HBI).
% nIter         = scalar integer, number of simulations (will be averaged
% to get 'recovered' parameter).
% nSub          = scalar integer, number of subjects.
%
% OUTPUTS:
% no outputs, just plots.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.
% Inspired by code from Algermissen, J. et al. 2024.
% ----------------------------------------------------------------------- %
%% 00a) Directories:
run("github_config.m") % Contains path to actual data folder; adjust with your onw path
dirs = [];
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.behav = fullfile(folderPath);
dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
dirs.sim     = fullfile(dirs.results, 'Simulations');

% Create new directories if necessary:
dirs.paramRecovery = fullfile(dirs.results, 'parameterRecovery');
if ~exist(dirs.paramRecovery, 'dir'); mkdir(dirs.paramRecovery); end

close all;
% ----------------------------------------------------------------------- %
%% 00b) Complete settings:

% Other settings:
selMod        = 23;
selMods = [14, 18, 22, 23, 33:40];
simType     = 'modSim'; % modSim realSim
parType     = 'hbi'; % lap hbi
nIter       = 100; % number simulations to be recovered
refit = input("Refit models? Y/N ", "s");
refit = strcmpi(refit, "y");
if ~refit
    refitIndividual = input("But refit individual models? Y/N ", "s");
    refitIndividual = strcmpi(refitIndividual, "y");
end
fprintf('Perform model recovery for M%02d, simType %s, parType %s \n', ...
    selMod, simType, parType);


% ----------------------------------------------------------------------- %
%% 00b) Add paths:
fprintf('Add paths\n')

addpath('/home/control/samnel/Documents/MATLAB/cbm-master/codes');
addpath(fullfile(dirs.root, 'Analyses/CBM_Scripts/models')); % models
addpath(fullfile(dirs.root, 'Analyses/CBM_Scripts')); % for sim_subj
addpath '/home/common/matlab/fieldtrip/qsub'; % qsubcellfun in Fieldtrip

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%% 00c) Prepare fitting simulated data:

% Priors
fprintf('Initialize priors\n')

priors{1} = struct('mean', [0 2], 'variance', [3 5]); % prior_model
priors{2} = struct('mean', [0 2 0], 'variance', [3 5 10]); % prior_model_goBias
priors{3} = struct('mean', [0 2 0 0], 'variance', [3 5 10 10]); % prior_model_fixedPavlov
priors{4} = struct('mean', [0 2 0 0], 'variance', [3 5 10 10]); % prior_model_dynamicPavlov
priors{5} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_fixedOmega
priors{6} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_dynamicOmega1
priors{7} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_fixedOmega
priors{8} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_dynamicOmega1
priors{9} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{10} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{11} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{12} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{13} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{14} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 10 3]); % prior_model_fixedPavlov
priors{15} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{16} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{17} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{18} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 10 3]); % prior_model_fixedPavlov
priors{19} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{20} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 10]); % prior_model_dynamicOmega2
priors{21} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{22} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{23} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 10 3]); % prior_model_fixedPavlov
priors{24} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{25} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{26} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{27} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{28} = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 10]); % prior_model_dynamicOmega2
priors{29} = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 10]); % prior_model_dynamicOmega2
priors{30} = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 10]); % prior_model_dynamicOmega2
priors{31} = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 10]); % prior_model_dynamicOmega2
priors{32} = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 3]); % prior_model_dynamicOmega2
priors{33} = struct('mean', [0 2 0 2 0 0], 'variance', [3 5 10 5 3 3]);
priors{34} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{35} = struct('mean', [0 2 0 2 0], 'variance', [3 5 10 5 3]); % prior_model_dynamicOmega2
priors{36} = struct('mean', [0 0 2 0 0], 'variance', [3 10 5 3 3]);
priors{37} = struct('mean', [0 0 0 3 0], 'variance', [3 10 3 5 3]);
priors{38} = struct('mean', [0 0 3 0], 'variance', [3 10 5 3]);
priors{39} = struct('mean', [2 0 0 0], 'variance', [5 10 10 3]); % prior_model_fixedPavlov
priors{40} = struct('mean', [2 0 0], 'variance', [5 10 10]); % prior_model_fixedPavlov


% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 01a) Load true parameters:

% Must have been generated with ChemControl_cbm_fit.

for i = 1:length(selMods)
    selMod = selMods(i);
fprintf('Load ground truth parameters for M%02d, parType %s\n', selMod, parType);

fprintf('Load parameters based on %s\n', parType);

if strcmp(parType, 'lap')
    lap_name_mod = fullfile(dirs.lap, sprintf('lap_mod%02d_allData.mat', selMod));
    fname = load(lap_name_mod);
    cbm = fname.cbm;
    trueParam = cbm.output.parameters;
elseif strcmp(parType, 'hbi')
    hbi_name_mod = fullfile(dirs.hbi, sprintf('hbi_mod_%02d_allData.mat', selMod));
    fname = load(hbi_name_mod);
    cbm = fname.cbm;
    trueParam = cbm.output.parameters{:};
end

% Initialize number parameters:
nParam = size(trueParam, 2);

% ----------------------------------------------------------------------- %
%% 01b) Load model simulations:

% Must have been generated with ChemControl_cbm_sim.

fprintf('Load modSim simulations for M%02d, parType %s\n', ...
    selMod, parType);

simFile     = fullfile(dirs.sim, sprintf('%s_mod%02d_%s_iter%04d.mat', ...
    simType, selMod, parType, nIter));

tmp         = load(simFile);
simulations = tmp.sim;
nSub        = size(simulations.stimuli, 1); % number subjects


%% 02a) Start Laplace approximation (LAP) via qsubeval jobs for each iteration:
refitNext = 0;
if ~refit && refitIndividual
    refitNext = input(sprintf("Refit M%02d? Y/N", selMod), "s");
    refitNext = strcmpi(refitNext, "y");
end


    if refit || refitNext
        for iIter = 1:nIter
        
            % ------------------------------------------------------------------- %
            % Prepare data for each subject for given iteration:
            simData = cell(nSub, 1);
            for iSub = 1:nSub
                simData{iSub}.stimuli = squeeze(simulations.stimuli(iSub, iIter, :, :));
                simData{iSub}.outcomes = squeeze(simulations.outcomes(iSub, iIter, :, :));
                simData{iSub}.actions = squeeze(simulations.actions(iSub, iIter, :, :));
            end
        
        % ------------------------------------------------------------------- %
            % Output name:
            fname_mod = fullfile(dirs.paramRecovery, sprintf('lap_mod%02d_iter%04d.mat', selMod, iIter));
            
            % ------------------------------------------------------------------ %
            % Fit with Laplace approximation sequentially
            % cbm_lap(simData, eval(sprintf('@ChemControl_mod%d', selMod)), priors{selMod}, fname_mod);
        
            % % -------------------------------------------------------------------
            % % Fit with Laplace approximation in parallel with qsubfeval on the
            % Donders HPC
            req_mem = 1024^3 * 5; % (bytes * KB * MB) * GB;
            req_etime = 60 * 60; % 30 minutes
            qsubfeval(@cbm_lap, simData, eval(sprintf('@ChemControl_mod%d', selMod)), priors{selMod}, fname_mod, 'memreq', req_mem, 'timreq', req_etime);
        end
    end
    % ------------------------------------------------------------------- %
    %% 02b) Read in LAP parameters fitted to simulated data:
    
    fprintf('Load fitted parameters for lap_mod%02d_iter%04d.mat \n', selMod, nIter);
    
    parType = 'lap';
    fittedParam = nan(nSub, nParam, nIter);
    
    for iIter = 1:nIter
        fname_lap = fullfile(dirs.paramRecovery, sprintf('lap_mod%02d_iter%04d.mat', selMod, iIter));
        
        % Poll for the file until it exists or until a maximum wait time is reached
        fileExists = false;
        
        while ~fileExists
            if exist(fname_lap, 'file')
                fileExists = true;
            else
                fprintf('File not found: %s. Retrying in 10 seconds...\n', fname_lap);
                pause(10); % Wait for 10 seconds before checking again
            end
        end
        
        % If file was found, proceed with loading
        if fileExists
            tmp = load(fname_lap);
            fittedParam(:, :, iIter) = tmp.cbm.output.parameters;
        else
            error('File not found after multiple attempts: %s', fname_lap);
        end
    end
    
    % ----------------------------------------------------------------------- %
    % ----------------------------------------------------------------------- %
    % ----------------------------------------------------------------------- %

%% 03a) Fit with Hierarchical Bayesian Inference (HBI) via qsubeval jobs for each iteration:

    if refit || refitNext
        % model code:
        models = {eval(sprintf('@ChemControl_mod%d', selMod))};
    
        for iIter = 1:nIter
            simData = cell(nSub, 1);
            for iSub = 1:nSub
                simData{iSub}.stimuli = squeeze(simulations.stimuli(iSub, iIter, :, :));
                simData{iSub}.outcomes = squeeze(simulations.outcomes(iSub, iIter, :, :));
                simData{iSub}.actions = squeeze(simulations.actions(iSub, iIter, :, :));
            end
            % ------------------------------------------------------------------- %
            % LAP fitted object:
            fcbm_maps = {fullfile(dirs.paramRecovery, sprintf('lap_mod%02d_iter%04d.mat', selMod, iIter))};
            % ------------------------------------------------------------------- %
            % HBI output file:
            fname_hbi = {fullfile(dirs.paramRecovery, sprintf('hbi_mod%02d_iter%04d.mat', selMod, iIter))};
            % ------------------------------------------------------------------- %
            % Fit with HBI in parallel with qsubfeval:
            req_mem = 1024^3 * 5;
            req_etime = 60 * 60; % 30 minutes to be safe
            qsubfeval(@cbm_hbi, simData, models, fcbm_maps, fname_hbi, 'memreq', req_mem, 'timreq', req_etime);
        end
    end

    % ----------------------------------------------------------------------- %
    %% 03b) Read in HBI parameters fitted to simulated data: 
    
    fprintf('Load fitted parameters for hbi_mod%02d_iterXXXX.mat \n', selMod);
    
    parType = 'hbi';
    fittedParam = nan(nSub, nParam, nIter);
    
    for iIter = 1:nIter
        fname_hbi = fullfile(dirs.paramRecovery, sprintf('hbi_mod%02d_iter%04d.mat', selMod, iIter));
        
        % Poll for the file until it exists or until a maximum wait time is reached
        fileExists = false;
        
        while ~fileExists
            if exist(fname_hbi, 'file')
                fileExists = true;
            else
                fprintf('File not found: %s. Retrying in 10 seconds...\n', fname_hbi);
                pause(10); % Wait for 10 seconds before checking again
            end
        end
        
        % If file was found, proceed with loading
        if fileExists
            tmp = load(fname_hbi);
            fittedParam(:, :, iIter) = tmp.cbm.output.parameters{:};
        else
            error('File not found after multiple attempts: %s', fname_hbi);
        end
    end

    
    % ----------------------------------------------------------------------- %
    % ----------------------------------------------------------------------- %
    % ----------------------------------------------------------------------- %
    % ----------------------------------------------------------------------- %

%% 04a) Descriptives and averaging:

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

param_names{7} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}'}; % DynamicOmega2Model
transform{7} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid'};

param_names{14} = {'\epsilon', '\rho', 'goBias', '\pi', '\alpha_{lr}'};
transform{14} = {'sigmoid', 'exp', '@(x) x', '@(x) x', 'sigmoid'};
param_names{18} = {'\epsilon', '\rho', 'goBias', '\pi', '\alpha_{lr}'}; % FixedPavlovModel
transform{18} = {'sigmoid', 'exp', '@(x) x', '@(x) x', 'sigmoid'};
param_names{21} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}'}; % DynamicOmega2Model
transform{21} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid'};
param_names{22} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{22} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid', 'sigmoid'};
param_names{23} = {'\epsilon', '\rho', 'goBias', '\pi', '\alpha_{lr}'}; % FixedPavlovModel
transform{23} = {'sigmoid', 'exp', '@(x) x', '@(x) x', 'sigmoid'};
param_names{24} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{24} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid', 'sigmoid'};
param_names{25} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{25} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid', 'sigmoid'};
param_names{26} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{26} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid', 'sigmoid'};
param_names{33} = {'\epsilon', '\rho', 'goBias','\beta_{\Omega}', 'thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{33} = {'sigmoid', 'exp', '@(x) x', 'exp', '@scaledSigmoid', 'sigmoid'};
param_names{34} = {'\epsilon', '\rho', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}'}; % DynamicOmega2Model
transform{34} = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid'};
param_names{35} = {'\epsilon', '\rho', 'goBias', '\beta_{\Omega}', 'thres_{\Omega}'}; % DynamicOmega2Model
transform{35} = {'sigmoid', 'exp', '@(x) x', 'exp', '@scaledSigmoid'};
param_names{36} = {'\epsilon', 'goBias','\beta_{\Omega}', 'thres_{\Omega}', '\alpha_{lr}'}; % DynamicOmega2Model
transform{36} = {'sigmoid', '@(x) x', 'exp', '@scaledSigmoid', 'sigmoid'};
param_names{37} = {'\epsilon', 'goBias', '\alpha_{\Omega}','\beta_{\Omega}', 'thres_{\Omega}'}; % DynamicOmega2Model
transform{37} = {'sigmoid', '@(x) x', 'sigmoid', 'exp', '@scaledSigmoid'};
param_names{38} = {'\epsilon', 'goBias', '\beta_{\Omega}', 'thres_{\Omega}'}; % DynamicOmega2Model
transform{38} = {'sigmoid', '@(x) x', 'exp', '@scaledSigmoid'};
param_names{39} = {'\rho', 'goBias', '\pi', '\alpha_{lr}'}; % FixedPavlovModel
transform{39} = {'exp', '@(x) x', '@(x) x', 'sigmoid'};
param_names{40} = {'\rho', 'goBias', '\pi'}; % FixedPavlovModel
transform{40} = {'exp', '@(x) x', '@(x) x'};

% Descriptives of parameters:
fittedParamAvgSub = squeeze(mean(fittedParam, 1));

% Average over iterations:
fittedParamAvgIter = squeeze(mean(fittedParam, 3));

% Descriptives of parameters:
fprintf('Mean parameter average over subjects: \n');
for iParam = 1:nParam
    fprintf('Parameter %s: M = %.02f, SD = %.02f, Min = %.02f, Max = %.02f\n', ...
        param_names{selMod}{iParam}, mean(fittedParamAvgSub(iParam, :)), std(fittedParamAvgSub(iParam, :)), ...
        min(fittedParamAvgSub(iParam, :)), max(fittedParamAvgSub(iParam, :)));
end

% ----------------------------------------------------------------------- %
%% 04b) Transform parameters:

fprintf('Transforming parameters for model %d \n', selMod);

for j = 1:length(transform{selMod})
    % Apply transformation to both trueParam and fittedParamAvgIter
    transformFunc = str2func(transform{selMod}{j});
    trueParam(:, nParam + j) = transformFunc(trueParam(:, j));
    fittedParamAvgIter(:, nParam + j) = transformFunc(fittedParamAvgIter(:, j));
end

% ----------------------------------------------------------------------- %
%% 04c) Detect outliers:

figure;
legendVec = {};

for iParam = 1:nParam
    % Plot true parameters
    subplot(1, nParam, iParam)
    hold on
    plot(trueParam(:, nParam + iParam), 'LineWidth', 2);
    plot(fittedParamAvgIter(:, nParam + iParam), 'LineWidth', 2);
    hold off
    legendVec{end+1} = sprintf('True Param %s', param_names{selMod}{iParam});
    legendVec{end+1} = sprintf('Fitted Param %s', param_names{selMod}{iParam});

    xlabel('Subject'); ylabel('Value');
    title(sprintf('True Param %s', param_names{selMod}{iParam}));
end

legend(legendVec);
sgtitle(sprintf("M%02d", selMod));
    
%trueParam(:, iParam);
%plot(trueParam(:, iParam));

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 05a) Select subjects:

fprintf('Select subjects\n');
fprintf('Select outlying subjects, and put them into a separate list')

invalidSubs = [];

fprintf('Exclude subjects %s \n', num2str(invalidSubs));
validSubs = setdiff(1:nSub, invalidSubs);

% ----------------------------------------------------------------------- %
%% 05b) Correlate true vs. fitted parameters:

for iParam = nParam+1:nParam*2
    if iParam <= nParam
        % Original parameters
        fprintf('Correlation parameter %s: Pearson r = %.03f, Spearman r = %.03f\n', ...
            param_names{selMod}{iParam}, ...
            corr(trueParam(validSubs, iParam), fittedParamAvgIter(validSubs, iParam), "type", "Pearson"), ...
            corr(trueParam(validSubs, iParam), fittedParamAvgIter(validSubs, iParam), "type", "Spearman"));
    else
        % Transformed parameters (after the original nParam)
        fprintf('Correlation parameter %s (transformed): Pearson r = %.03f, Spearman r = %.03f\n', ...
            param_names{selMod}{iParam-nParam}, ...  % Reference the original parameter name for the transformed version
            corr(trueParam(validSubs, iParam), fittedParamAvgIter(validSubs, iParam), "type", "Pearson"), ...
            corr(trueParam(validSubs, iParam), fittedParamAvgIter(validSubs, iParam), "type", "Spearman"));
    end
end

% Save source data files:
fullFileName = fullfile(dirs.root, 'Log/Plots/trueParam.csv');
csvwrite(fullFileName, trueParam(validSubs, :));
fullFileName = fullfile(dirs.root, 'Log/Plots/fittedParam.csv');
csvwrite(fullFileName, fittedParamAvgIter(validSubs, :));



%% 05c) Correlate between parameters
corrMatrix = zeros(nParam, nParam);
param_labels = param_names{selMod};

for iParam1 = 1:nParam
    for iParam2 = 1:nParam 
        
        x = fittedParamAvgIter(:, nParam + iParam1);
        y = fittedParamAvgIter(:, nParam + iParam2);
        
        corrMatrix(iParam1, iParam2) = round(corr(x, y), 2);
    end
end
figure;
% Create the heatmap
h = heatmap(corrMatrix, 'XData', param_labels, 'YData', param_labels, ...
            'XLabel', 'Parameter', 'YLabel', 'Parameter', ...
            'Title', sprintf('Correlation Matrix (Model %02d)', selMod), ...
            'ColorLimits', [-1 1], 'Colormap', parula);

% Move the x-tick labels to the top
ax = struct(h);  % Access the underlying axes
ax.Axes.XAxisLocation = 'top';  % Move x-axis ticks to the top
end
% ----------------------------------------------------------------------- %
%% 05c) Plot true vs. fitted parameters:

% Choose parameter:

% for iParam = 1:nParam * 2
%     % Plotting settings:
%     fontSize        = 36;
%     lineWidth       = 4;
%     markerSize      = 35; % 25 
% 
%     % Select variables:
%     x = trueParam(:, iParam);
%     y = fittedParamAvgIter(:, iParam);
% 
%     % Start plot:
%     close all
%     figure('Position', [100, 100 1200 1200]);
% 
%     % Plot scatterplot
%     plot(x, y, 'k.', 'MarkerSize', markerSize)
% 
%     hold on;
%     % Add regression line:
%     b = polyfit(x(validSubs), y(validSubs), 1); % fit regression
%     Yhat = polyval(b, x(validSubs)); % predict
%     plot(x(validSubs), Yhat, 'r-', 'Linewidth', lineWidth)
% 
%     % Add identity line:
%     plot([-100 100], [-100 100], 'k--', 'Linewidth', 1)
% 
%     % Axis limits:
%     minLim = min([x; y]);
%     maxLim = max([x; y]);
%     axisLim = [minLim - 0.1*maxLim maxLim + 0.1*maxLim];
%     xlim(axisLim);
%     ylim(axisLim);
% 
%     % Settings:
%     set(gca, 'Fontsize', fontSize, 'Linewidth', lineWidth);
%     if iParam < 5
%         title(sprintf('Untransformed Parameter %s', param_names{selMod}{iParam}));
%     else
%         title(sprintf('Transformed Parameter %s', param_names{selMod}{-4 + iParam}));
%     end
% 
%     xlabel('True parameter');
%     ylabel('Fitted parameter');
% 
%     % Save:
%     figName = fullfile(dirs.root, sprintf('Log/Plots/paramRecovery_%s_mod%02d_param%02d', ...
%         parType, selMod, iParam));
%     if ~isempty(invalidSubs); figName = sprintf('%s_without%s', figName, strjoin(string(invalidSubs), '_')); end
%     saveas(gcf, [figName '.png']);
% 
%     % Close:
%     pause(1)
%     close gcf
% end