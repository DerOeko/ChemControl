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

% priors{1} = struct('mean', [0 2], 'variance', [3 5]); % prior_model
% priors{2} = struct('mean', [0 2 0], 'variance', [3 5 10]); % prior_model_goBias
% priors{3} = struct('mean', [0 2 0 0], 'variance', [3 5 10 10]); % prior_model_fixedPavlov
% priors{4} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_dynamicPavlov
% priors{5} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{6} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{7} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{8} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{9} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{10} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{11} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{12} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{13} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{14} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{15} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{16} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{17} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{18} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{19} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{20} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{21} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{22} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{23} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{24} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{25} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{26} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{27} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{28} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{29} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{30} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{31} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{32} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{33} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{34} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{35} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{36} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{37} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{38} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{39} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{40} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{41} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{42} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{43} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{44} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{45} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{46} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{47} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{48} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{49} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{50} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{51} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{52} = struct('mean', [0 2 0 2 0 0], 'variance', [3 5 10 5 3 3]);
% priors{53} = struct('mean', [0 2 0 0 2 0 ], 'variance', [3 5 10 3 5 3]);
% priors{54} = struct('mean', [0 2 0 2 0 ], 'variance', [3 5 10 5 3]);
% priors{55} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]); % prior_model_dynamicPavlov
% priors{56} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{57} = struct('mean', [0 2 0 0 2 0 0], 'variance', [3 5 10 3 5 3 3]);
% priors{58} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{59} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{60} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{61} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{62} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{63} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{64} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{65} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{66} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{67} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{68} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{69} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{70} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{71} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{72} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{73} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{74} = struct('mean', [0 2 0 0 2], 'variance', [3 5 10 3 5]);
% priors{75} = struct('mean', [0 2 0 0 2 2], 'variance', [3 5 10 3 5 5]);
% priors{76} = struct('mean', [0 2 0 0 2 2 0], 'variance', [3 5 10 3 5 5 10]);
% priors{77} = struct('mean', [0 2 0 0 2], 'variance', [3 5 10 3 5]);
% priors{78} = struct('mean', [0 2 0 0 2 2], 'variance', [3 5 10 3 5 5]);
% priors{79} = struct('mean', [0 2 0 0 2 2 0], 'variance', [3 5 10 3 5 5 10]);
% priors{80} = struct('mean', [0 0 0 2 0], 'variance', [3 10 3 5 3]);
% priors{81} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{82} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{83} = struct('mean', [0 0 2 0 0 2 0], 'variance', [3 3 5 10 3 5 3]);
% priors{84} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{85} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{86} = struct('mean', [0 2 0 0 2 0], 'variance', [3 5 10 3 5 3]);
% priors{87} = struct('mean', [2 0 2 0 2 0], 'variance', [5 3 5 3 5 3]);
priors{1} = struct('mean', [0 2], 'variance', [3 5]);
priors{2} = struct('mean', [0 2 0], 'variance', [3 5 10]);
priors{3} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]);
priors{4} = struct('mean', [0 2 0 0], 'variance', [3 5 10 3]);
priors{5} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 3 3]);
priors{6} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 3 3]);
priors{7} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 3 3]);
priors{8} = struct('mean', [0 2 0 0 0], 'variance', [3 5 10 3 3]);

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
    parameters = randn(1, 30);
    F1 = eval(sprintf('ChemControl_mod%d(parameters, subj1)', iMod));
    fprintf('Model %02d: loglik = %f\n', iMod, F1);
end

fprintf(">>>>>  Test with extreme values\n")
% b) Extreme parameter values:
for iMod = selMods
    parameters = ones(30) * 10;
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
control_types = {'allData'};

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
    for iMod = 1:length(selMods)
        selMod = selMods(iMod);

        outputFileName = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', selMod, ctype));
        % Check if the model has already been fitted:
        if exist(outputFileName, 'file')
            if refitAll
                fprintf('Refitting model %02d with Laplace approximation for %s data\n', selMod, ctype);
                cbm_lap(d, eval(sprintf('@ChemControl_mod%d', selMod)), priors{selMod}, outputFileName);
            else
                fprintf('Model %02d already fitted for %s data. Refit? Y/N [N]: ', selMod, ctype);
                choice = input('', 's');
                if isempty(choice)
                    choice = 'N';
                end
                if upper(choice) == 'Y'
                    fprintf('Refitting model %02d with Laplace approximation for %s data\n', selMod, ctype);
                    cbm_lap(d, eval(sprintf('@ChemControl_mod%d', selMod)), priors{selMod}, outputFileName);
                else
                    fprintf('Skipping refit for model %02d (%s)\n', selMod, ctype);
                end
            end
        else
            fprintf('Fit model %02d with Laplace approximation for %s data\n', selMod, ctype);
            cbm_lap(d, eval(sprintf('@ChemControl_mod%d', selMod)), priors{selMod}, outputFileName);
        end
    end

    % Step 2: Evaluate LAP Fit
    fprintf('Evaluating LAP fit for %s data\n', ctype);
    % Initialize output file names
    fname_mod = cell(length(selMods), 1);
    logModEvi = nan(nSub, length(selMods));
    for iMod = 1:length(selMods)
        fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d_%s.mat', selMods(iMod), ctype));
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
    for iMod = 1:length(selMods)
        plot(1:nSub, logModEvi(:, iMod), 'LineWidth', 2, 'DisplayName', sprintf('Model %02d', selMods(iMod)));
        legendVec = [legendVec; sprintf('Model %02d', selMods(iMod))];
    end
    xlabel('Subjects'); ylabel('Log model evidence');
    legend(legendVec);

    % t-test between 2 models
    for iMod1 = 1:length(selMods)
        for iMod2 = 1:length(selMods)
            selMod1 = selMods(iMod1);
            selMod2 = selMods(iMod2);
            if iMod1 ~= iMod2
                fprintf('Test log model-evidence of models %d and %d against each other for %s: \n', selMod1, selMod2, ctype);
                [~, p, ~, STATS] = ttest(logModEvi(:, iMod1), logModEvi(:, iMod2));
                fprintf('t(%d) = %.3f, p = %.03f\n', STATS.df, STATS.tstat, p);
            end
        end
    end

    % Model frequency
    modRange = 1:length(selMods);
    modFreq = nan(nSub, 1);
    for iSub = 1:nSub
        [~, I] = max(logModEvi(iSub, modRange));
        modFreq(iSub) = I;
    end

    % Model frequency per model
    fprintf('Model frequency per model for %s: \n', ctype);
    for iMod = modRange
        selMod = selMods(iMod);
        idx = find(modFreq == iMod);
        fprintf('Model M%02d is best for %02d subjects: %s\n', selMod, length(idx), mat2str(idx));
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
    for iMod = 1:length(selMods)
        selMod = selMods(iMod);
        idx = find(subResp == iMod);
        fprintf('Model M%02d is responsible for %02d subjects: %s\n', selMod, length(idx), mat2str(idx));
    end

    % Exceedance probability:
    fprintf('Exceedance probability for %s:\n', ctype);
    disp(cbm.output.exceedance_prob)

    % Protected exceedance probability:
    fprintf('Protected exceedance probability for %s:\n', ctype);
    disp(cbm.output.protected_exceedance_prob)
end

%% AIC, BIC, Model Likelihood comparisons

logs = -logModEvi;
AICvalues = zeros(nSub, length(selMods ));
BICvalues = zeros(nSub, length(selMods ));
MDLvalues = zeros(nSub, length(selMods));

for iMod = 1:length(selMods)
    nParam = size(cbm.output.parameters{iMod}, 2);

    for iSub = 1:nSub
        L = logModEvi(iSub, iMod);
        n = 360;
        AICvalues(iSub, iMod) = 2*nParam - 2 * L;
        BICvalues(iSub, iMod) = nParam*log(n) - 2 * L;
        MDLvalues(iSub, iMod) = mdl_rl(nParam, n, -L);
    end
end
      

% Define model names
modelNames = cell(length(selMods), 1); % Update with your actual model names
for iMod = 1:length(selMods)
    selMod = selMods(iMod);
    modelNames{iMod} = sprintf('M%02d', selMod);
end

x = 1:length(selMods);

figure;
hold on;
bar(x, mean(AICvalues, 1), 0.5, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5)

for iSub = 1:nSub
    plot(x, AICvalues(iSub, :), '-o', 'LineWidth', 1);
end
xlabel('Model');
ylabel('AIC Value');
title('AIC Values Across Models for Each Subject');
xticks(x);
xticklabels(modelNames);
grid on;

% Optionally, add a legend if the number of subjects is small
% legend(arrayfun(@(s) sprintf('Subject %d', s), 1:nSub, 'UniformOutput', false), 'Location', 'BestOutside');

hold off;

figure;
boxplot(AICvalues, 'Labels', modelNames);
xlabel('Model');
ylabel('AIC Value');
title('AIC Values Distribution Across Subjects');
grid on;

figure;
hold on;
bar(x, mean(BICvalues, 1), 0.5, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5)

for iSub = 1:nSub
    plot(x, BICvalues(iSub, :), '-o', 'LineWidth', 1);
end

xlabel('Model');
ylabel('BIC Value');
title('BIC Values Across Models for Each Subject');
xticks(x);
xticklabels(modelNames);
grid on;

% Optionally, add a legend if the number of subjects is small
% legend(arrayfun(@(s) sprintf('Subject %d', s), 1:nSub, 'UniformOutput', false), 'Location', 'BestOutside');

hold off;


figure;
hold on;
bar(x, mean(MDLvalues, 1), 0.5, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5)

for iSub = 1:nSub
    plot(x, MDLvalues(iSub, :), '-o', 'LineWidth', 1);
end
xlabel('Model');
ylabel('MDL Value');
title('MDL Values Across Models for Each Subject');
xticks(x);
xticklabels(modelNames);
grid on;

% Optionally, add a legend if the number of subjects is small
% legend(arrayfun(@(s) sprintf('Subject %d', s), 1:nSub, 'UniformOutput', false), 'Location', 'BestOutside');

hold off;

figure;
boxplot(BICvalues, 'Labels', modelNames);
xlabel('Model');
ylabel('BIC Value');
title('BIC Values Distribution Across Subjects');
grid on;

sumAIC = sum(AICvalues, 1);
sumBIC = sum(BICvalues, 1);
% Plot the bar chart
figure;
hold on;
% Create grouped bar chart
barData = [sumAIC; sumBIC]';
b = bar(barData, 'grouped');

% Set colors for AIC and BIC bars
b(1).FaceColor = [0.2, 0.6, 0.5]; % Teal for AIC
b(2).FaceColor = [0.8, 0.4, 0.4]; % Reddish for BIC

% Add labels and title
xlabel('Models');
ylabel('Sum of AIC/BIC Values');
title('Comparison of AIC and BIC Across Models');

% Set x-axis tick labels to model names
xticks(1:length(modelNames));
xticklabels(modelNames);

% Add legend
legend({'AIC', 'BIC'}, 'Location', 'Best');

% Improve appearance with gridlines
grid on;

% Display values on top of bars
for k = 1:length(b)
    % Get X and Y data of bars for placing text
    xData = b(k).XEndPoints;
    yData = b(k).YEndPoints;
    labels = string(round(b(k).YData, 1)); % Round the values for display
    text(xData, yData, labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10);
end

hold off;
% END OF FILE.

% modVec = selMods;
% fname_hbi = fullfile(dirs.hbi, sprintf('hbi_mod%s_%s.mat', num2str(modVec, '_%02d'), "allData"));
% f_hbi = load(fname_hbi);
% cbm = f_hbi.cbm;
% nSub = size(cbm.output.parameters{1}, 1);
% responsibility = cbm.output.responsibility;
% 
% subResp = nan(nSub, 1);
% for iSub = 1:nSub
%     subMax = max(responsibility(iSub, :));
%     subResp(iSub) = find(responsibility(iSub, :) == subMax);
% end
% tabulate(subResp);

    