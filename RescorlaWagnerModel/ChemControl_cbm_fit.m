% ChemControl_fit.m

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
selMods = 1:nMod;
selMods = [1:5 11 12 16 22];
fprintf('Fit %d models\n', length(selMods));

% ----------------------------------------------------------------------- %
%% 00c) Paths:
fprintf('Add paths\n');

% Add paths: 
addpath('/home/control/samnel/Documents/MATLAB/cbm-master/codes'); % CBM toolbox
addpath(fullfile(dirs.root, 'behavioral_study', 'scripts', 'matlab_scripts', 'RescorlaWagnerModel', 'models')); % models

% ----------------------------------------------------------------------- %
%% 00d) Load and extract data, priors, output name:

fprintf('Load data\n');
inputFile = fullfile(dirs.results, 'ChemControl_cbm_inputData.mat');
fdata = load(inputFile);
allData_d = fdata.data;
nSub = length(allData_d);

%% 01a) Alternative a) Fitting with Stan
%% 01b) Alternative b) Fitting with CBM
% Priors:
fprintf('Initialize priors\n')

priors{1} = struct('mean', [0 2], 'variance', [3 5]); % prior_model
priors{2} = struct('mean', [0 2 0], 'variance', [3 5 10]); % prior_model_goBias
priors{3} = struct('mean', [0 2 0 0], 'variance', [3 5 10 10]); % prior_model_fixedPavlov
priors{4} = struct('mean', [0 2 0 0], 'variance', [3 5 10 10]); % prior_model_dynamicPavlov
priors{5} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_fixedOmega
priors{6} = struct('mean', [0 2 0 0 0 2], 'variance', [3 5 10 3 3 5]); % prior_model_dynamicOmega1
priors{7} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{8} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega3
priors{9} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega4
priors{10} = struct('mean', [0 2 0 2 0 0], 'variance', [3 5 10 5 3 10]); % prior_model_dynamicOmega5
priors{11} = struct('mean', [0 2 0], 'variance', [3 5 10]); % prior_model_dynamicOmega6
priors{12} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 3 3]); % prior_model_dynamicOmega5
priors{13} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{14} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{15} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{16} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]); % prior_model_dynamicOmega2
priors{17} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{18} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{19} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{20} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{21} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]); % prior_model_dynamicOmega2
priors{22} = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 3]); % prior_model_dynamicOmega2

% Output names:
fprintf("Initialize output file names\n")

fname_mod = cell(length(selMods), 1);
for iMod = selMods
    fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d.mat', iMod));
end

% ----------------------------------------------------------------------- %
% 01b) Check models in dry run:

fprintf('Test models (dry run)\n')
subj1 = allData_d{1};

fprintf(">>>>>  Test with random values\n")
% a) Random parameter values:
for iMod = selMods
    parameters = randn(1, 8);
    F1 = eval(sprintf('ChemControl_mod%d(parameters, subj1)', iMod));
    fprintf('Model %02d: loglik = %f\n', iMod, F1);
end

fprintf(">>>>>  Test with extreme values\n")
% b) Extreme parameter values:
for iMod = selMods
    parameters = [-10 10 -10 10 10 10 10 10];
    F1 = eval(sprintf('ChemControl_mod%d(parameters, subj1)', iMod));
    fprintf('Model %02d: loglik = %f\n', iMod, F1);
end


%% 2 Fitting to subsets of data separately
% 2.1 extracting high and low control data
hc_d = cell(1, nSub);
lc_d = cell(1, nSub);
nBlocks = size(allData_d{1}.controllability, 1);

for iSub = 1:nSub
    d = allData_d{iSub};
    % Identify high control (hc) blocks
    hc_idx = find(d.controllability(:, 1) == 1);
    
    % Identify low control (lc) blocks where controllability is 0 and isYoked is 0
    lc_idx = find(d.controllability(:, 1) == 0);

    % Extract the data for the high control blocks and store it
    hc_d{iSub} = structfun(@(x) x(hc_idx, :), d, 'UniformOutput', false);
    lc_d{iSub} = structfun(@(x) x(lc_idx, :), d, 'UniformOutput', false);
end

%% 2.2 Fit with separate datasets and control types

% Define control types
control_types = {'allData', 'hc', 'lc'};

% Ask if all models should be refit:
refitAllQuery = 'Do you want to refit all models? Y/N [N]: ';
refitAllChoice = input(refitAllQuery, 's');
if isempty(refitAllChoice)
    refitAllChoice = 'N';
end
refitAll = upper(refitAllChoice) == 'Y';

% Unified processing loop for all control types
for ctype = control_types
    ctype = ctype{1}; % Extract string from cell
    fprintf('Processing %s data\n', ctype);
    
    % Load the appropriate data
    d = eval(sprintf('%s_d', ctype));

    % Step 1: Fit Laplace Approximation
    fprintf('Fitting Laplace Approximation for %s data\n', ctype);
    for iMod = selMods
        outputFileName = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, ctype));
        
        % Check if the model has already been fitted:
        if exist(outputFileName, 'file')
            if refitAll
                fprintf('Refitting model %02d with Laplace approximation for %s data\n', iMod, ctype);
                cbm_lap(d, eval(sprintf('@ChemControl_mod%d', iMod)), priors{iMod}, outputFileName);
            else
                fprintf('Model %02d already fitted for %s data. Refit? Y/N [N]: ', iMod, ctype);
                choice = input('', 's');
                if isempty(choice)
                    choice = 'N';
                end
                if upper(choice) == 'Y'
                    fprintf('Refitting model %02d with Laplace approximation for %s data\n', iMod, ctype);
                    cbm_lap(d, eval(sprintf('@ChemControl_mod%d', iMod)), priors{iMod}, outputFileName);
                else
                    fprintf('Skipping refit for model %02d (%s)\n', iMod, ctype);
                end
            end
        else
            fprintf('Fit model %02d with Laplace approximation for %s data\n', iMod, ctype);
            cbm_lap(d, eval(sprintf('@ChemControl_mod%d', iMod)), priors{iMod}, outputFileName);
        end
    end

    % Step 2: Evaluate LAP Fit
    fprintf('Evaluating LAP fit for %s data\n', ctype);
    % Initialize output file names
    fname_mod = cell(length(selMods), 1);
    logModEvi = nan(nSub, length(selMods));
    for iMod = selMods
        fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, ctype));
        fname = load(fname_mod{iMod});
        cbm = fname.cbm;
        logModEvi(:, iMod) = cbm.output.log_evidence;
    end

    % Mean per model
    fprintf('Mean log-model evidence per model for %s:\n', ctype);
    disp(round(mean(logModEvi, 1)))

    % Sum per model
    fprintf('Summed log-model evidence per model for %s:\n', ctype);
    disp(round(sum(logModEvi, 1)))

    % Plot log model-evidence per subject
    figure('unit', 'normalized', 'outerposition', [0 0 1 1]); hold on
    legendVec = [];
    for iMod = selMods
        plot(1:nSub, logModEvi(:, iMod), 'LineWidth', 2, 'DisplayName', sprintf('Model %02d', iMod));
        legendVec = [legendVec; sprintf('Model %02d', iMod)];
    end
    xlabel('Subjects'); ylabel('Log model evidence');
    legend(legendVec);

    % t-test between 2 models
    for iMod1 = selMods
        for iMod2 = selMods
            if iMod1 ~= iMod2
                fprintf('Test log model-evidence of models %d and %d against each other for %s: \n', iMod1, iMod2, ctype);
                [~, p, ~, STATS] = ttest(logModEvi(:, iMod1), logModEvi(:, iMod2));
                fprintf('t(%d) = %.3f, p = %.03f\n', STATS.df, STATS.tstat, p);
            end
        end
    end

    % Model frequency
    modRange = selMods;
    modFreq = nan(nSub, 1);
    for iSub = 1:nSub
        [~, I] = max(logModEvi(iSub, modRange));
        modFreq(iSub) = I;
    end

    % Model frequency per model
    fprintf('Model frequency per model for %s: \n', ctype);
    for iMod = modRange
        idx = find(modFreq == iMod);
        fprintf('Model M%02d is best for %02d subjects: %s\n', iMod, length(idx), mat2str(idx));
    end

    % Step 3: Fit Hierarchical Bayesian Inference (HBI) for Separate Models
    fprintf('Fitting HBI for %s data\n', ctype);
    for iMod = selMods
        fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod_%02d_%s.mat', iMod, ctype));
        
        % Check if the HBI model has already been fitted:
        if exist(fname_hbi, 'file')
            if refitAll
                fprintf('Refitting HBI for model %02d (%s)\n', iMod, ctype);
                % Format data, models, lap output files, and fit
                models = {eval(sprintf('@ChemControl_mod%d', iMod))};
                fcbm_maps = {fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, ctype))};
                cbm_hbi(d, models, fcbm_maps, {fname_hbi});
            else
                fprintf('HBI for model %02d already fitted for %s. Refit? Y/N [N]: ', iMod, ctype);
                choice = input('', 's');
                if isempty(choice)
                    choice = 'N';
                end
                if upper(choice) == 'Y'
                    fprintf('Refitting HBI for model %02d (%s)\n', iMod, ctype);
                    % Format data, models, lap output files, and fit
                    models = {eval(sprintf('@ChemControl_mod%d', iMod))};
                    fcbm_maps = {fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, ctype))};
                    cbm_hbi(d, models, fcbm_maps, {fname_hbi});
                else
                    fprintf('Skipping refit for HBI model %02d (%s)\n', iMod, ctype);
                end
            end
        else
            fprintf('Fit HBI for model %02d (%s)\n', iMod, ctype);
            % Format data, models, lap output files, and fit
            models = {eval(sprintf('@ChemControl_mod%d', iMod))};
            fcbm_maps = {fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', iMod, ctype))};
            cbm_hbi(d, models, fcbm_maps, {fname_hbi});
        end
    end

    % Step 4: Prepare Hierarchical Bayesian inference (cbm_hbi) across models for each control type
    fprintf('Prepare HBI for comparing models for %s\n', ctype);

    % Create names of input files:
    models = cell(length(selMods), 1);
    fcbm_maps = cell(length(selMods), 1);
    for iMod = 1:length(selMods)
        selMod = selMods(iMod);
        models{iMod} = eval(sprintf('@ChemControl_mod%d', selMod));
        fcbm_maps{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', selMod, ctype));
    end

    % File address for saving the output
    fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s_%s', num2str(selMods, '_%02d'), ctype));
    fprintf('Output file will be %s\n', fname_hbi);
    fname_hbi = [fname_hbi, '.mat'];
    
    % Check if the HBI model across models has already been fitted:
    if exist(fname_hbi, 'file')
        if refitAll
            fprintf('Refitting HBI for comparing models for %s\n', ctype);
            cbm_hbi(d, models, fcbm_maps, fname_hbi);
        else
            fprintf('HBI for comparing models already fitted for %s. Refit? Y/N [N]: ', ctype);
            choice = input('', 's');
            if isempty(choice)
                choice = 'N';
            end
            if upper(choice) == 'Y'
                fprintf('Refitting HBI for comparing models for %s\n', ctype);
                cbm_hbi(d, models, fcbm_maps, fname_hbi);
            else
                fprintf('Skipping refit for HBI comparing models for %s\n', ctype);
            end
        end
    else
        fprintf('Fit models %s with HBI for %s\n', num2str(selMods, '_%02d'), ctype);
        cbm_hbi(d, models, fcbm_maps, fname_hbi);
        fprintf("Finished fitting for %s :]\n", ctype);
    end

    % Additionally run protected exceedance probability (including null model)
    fprintf('Re-run models %s with HBI including null model for %s\n', num2str(selMods, '_%02d'), ctype);
    
    % Check if the HBI null model has already been fitted:
    if exist(fname_hbi, 'file')
        if refitAll
            fprintf('Refitting HBI null model for %s\n', ctype);
            cbm_hbi_null(d, fname_hbi);
        else
            fprintf('HBI null model already fitted for %s. Refit? Y/N [N]: ', ctype);
            choice = input('', 's');
            if isempty(choice)
                choice = 'N';
            end
            if upper(choice) == 'Y'
                fprintf('Refitting HBI null model for %s\n', ctype);
                cbm_hbi_null(d, fname_hbi);
            else
                fprintf('Skipping refit for HBI null model for %s\n', ctype);
            end
        end
    else
        fprintf('Fit HBI null model for %s\n', ctype);
        cbm_hbi_null(d, fname_hbi);
        fprintf('Finished fitting including null model for %s! :]\n', ctype);
    end

    % Evaluate HBI fit
    fprintf('Evaluating HBI fit for %s data\n', ctype);
    
    % Create output name
    modVec = selMods;
    fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s_%s.mat', num2str(modVec, '_%02d'), ctype));

    % Load model:
    f_hbi = load(fname_hbi);
    cbm = f_hbi.cbm;
    nSub = size(cbm.output.parameters{1}, 1);

    % Model frequency:
    fprintf('Model frequency for %s\n', ctype);
    disp(cbm.output.model_frequency)

    % Responsibility per subject:
    fprintf('Maximal model responsibility for %s\n', ctype);
    responsibility = cbm.output.responsibility;

    % Per subject:
    subResp = nan(nSub, 1);
    for iSub = 1:nSub
        subMax = max(responsibility(iSub, :));
        subResp(iSub) = find(responsibility(iSub, :) == subMax);
    end
    tabulate(subResp);

    % Per model:
    for iMod = selMods
        idx = find(subResp == iMod);
        fprintf('Model M%02d is responsible for %02d subjects: %s\n', iMod, length(idx), mat2str(idx));
    end

    % Exceedance probability:
    fprintf('Exceedance probability for %s:\n', ctype);
    disp(cbm.output.exceedance_prob)

    % Protected exceedance probability:
    fprintf('Protected exceedance probability for %s:\n', ctype);
    disp(cbm.output.protected_exceedance_prob)
end
% END OF FILE.

modVec = selMods;
fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s_%s.mat', num2str(modVec, '_%02d'), "allData"));
f_hbi = load(fname_hbi);
cbm = f_hbi.cbm;
nSub = size(cbm.output.parameters{1}, 1);
responsibility = cbm.output.responsibility;

subResp = nan(nSub, 1);
for iSub = 1:nSub
    subMax = max(responsibility(iSub, :));
    subResp(iSub) = find(responsibility(iSub, :) == subMax);
end
tabulate(subResp);

    