% ChemControl_testingModels.m

% This is an interactive script---execute it step-by-step.
% It tests a series of models and averages their runs.
% Made for purposes of seeing whether the retrieved parameters return good
% estimates.
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

% Directories:
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.results = fullfile(dirs.root, 'Log', 'Behavior', 'Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
dirs.models         = fullfile(dirs.root, 'modelSimulations');
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

% Input file
inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
data = load(inputFile).data;
fprintf('Loaded inputFile with %d subjects\n', size(data, 2));

nMod = length(dir(fullfile(dirs.models, "*.m")));
dataType = 'allData';

%% 01) LOAD AND SAVE:
%% 01a) Alternative A: Load model fitted with LAP:

parType = 'lap';
selMod = 1;
selMods = [3, 80, 81, 82, 83, 84, 85, 86];

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
%% 01b) Alternative B: Load model fitted with HBI:

parType = 'hbi';
% extract
for iMod = 1:length(selMods)
    idx = selMods(iMod); 
    % load
    fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d_%s.mat', idx, dataType));
    fprintf("Load model %d fit with HBI\n", idx);
    load(fname_hbi)
    fprintf("Extract model %d fit with HBI\n", idx);
    groupParams{iMod} = cbm.output.group_mean{:};
end

%% 01c) Check whether modSim and mod produce same loglikehood for same data.
% for iMod = 1:nMod
%     subj = sim_subj(8, 40);  % Simulate subject data
%     % Execute model simulation and retrieve output
%     out = eval(sprintf("ChemControl_mod%d_modSim(groupParams{iMod}, subj)", iMod));
%     loglik1 = out.loglik;  % Loglikelihood from simulation
% 
%     % Assign simulated actions back to subject for model re-evaluation
%     subj.actions = out.actions;
%     subj.outcomes = out.outcomes;
% 
%     % Calculate loglikelihood with the standard model function
%     loglik2 = eval(sprintf("ChemControl_mod%d(groupParams{iMod}, subj)", iMod));
% 
%     % Check if the loglikelihood values are the same
%     if loglik1 == loglik2
%         fprintf("Model %02d is producing the same loglikelihood for the same data and parameters.\n", iMod);
%     else
%         fprintf("Model %02d is not producing the same loglikelihood! Expected %f, got %f.\n", iMod, loglik1, loglik2);
%     end
% end


%% 02) SIMULATE

%% 02a) Run settings
%% Schedules, runs, trials, blocks, states, and controllability settings
nRuns = nSub*2;
nTrials = 40;
nBlocks = 8;
nStates = 4;
nSchedules = 11;
controllabilitySchedules = [
    1, 2, 2, 1, 2, 1, 1, 2;
    2, 1, 2, 2, 1, 2, 1, 1;
    1, 2, 1, 2, 2, 1, 2, 1;
    1, 1, 2, 1, 2, 2, 1, 2;
    2, 1, 1, 2, 1, 2, 2, 1;
    1, 2, 1, 1, 2, 1, 2, 2;
    2, 1, 2, 1, 1, 2, 1, 2;
    2, 2, 1, 2, 1, 1, 2, 1;
    1, 1, 2, 2, 1, 2, 1, 2;
    2, 2, 1, 1, 2, 1, 2, 1;
    1, 2, 1, 2, 1, 2, 1, 2;
];

cPs = zeros(nSchedules, nTrials*nBlocks);
for iSchedule = 1:nSchedules
    for iB = 1:nBlocks
        index = mod(iB -1, 8)+1;
        if controllabilitySchedules(iSchedule, index) == 1
            cPs(iSchedule, ((iB-1)*nTrials)+1:iB*nTrials+1)= 0.8;
        else
            cPs(iSchedule, ((iB-1)*nTrials)+1:iB*nTrials+1)=0.2;
        end
    end
end

fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig1)
sgtitle(sprintf("Learning curves for different models in high control trials for %i runs, %i blocks, %i trials", nRuns, nBlocks, nTrials))

fig2 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig2)
sgtitle(sprintf("Learning curves for different models in low control trials for %i runs, %i blocks, %i trials", nRuns, nBlocks, nTrials))

fig3 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig3)
sgtitle(sprintf("Learning curves for different models in yoked low control trials for %i runs, %i blocks, %i trials", nRuns, nBlocks, nTrials))

fig4 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig4)
sgtitle(sprintf("Proportion of staying with the same action when encountering the same stimulus as a function of previous action and previous received outcome"))

% Preallocate cell arrays for all outputs
allHCcell = cell(nRuns, 1);
allLCcell = cell(nRuns, 1);
allYCcell = cell(nRuns, 1);

% Preallocate cell arrays for all actions
allActions = cell(nRuns, 1);

% Preallocate cell arrays for all outcomes
allOutcomes = cell(nRuns, 1);

% Preallocate cell arrays for all stimuli
allStimuli = cell(nRuns, 1);

omegas = {};
omegasStim = zeros(nTrials * nBlocks, 4);
qs = {};
svs = {};
allOmegaStim = {};
i = 0;
j= 0;
k = 0;
l = 0;
schedule_counter = 0;
% 02b) Run simulation and plot
for iMod = 1:length(selMods)
    idx = selMods(iMod);
    parameters = groupParams{iMod};
    i = i+1;
    % Test
    subj = sim_subj(nBlocks, nTrials);
    out = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", idx));
    if isfield(out, 'omegas') 
        j= j + 1;
        modelsWithOmegas(j) = sprintf("M%02d", idx);
    end
    
    if isfield(out, 'omegasStim') 
        k= k + 1;
        modelsWithOmegasStim(k) = sprintf("M%02d", idx);
    end
    if isfield(out, 'qs')
        l = l + 1;
        modelsWithQs(l) = sprintf("M%02d", idx);
    end

    for iRun = 1:nRuns
        subj = sim_subj(nBlocks, nTrials);
        out = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", idx));
        % Collect omegas if they are a field in the output
        if isfield(out, 'omegas') && subj.selected_schedule_idx == 1
            % if size(out.omegas, 3) == 4
            %     reshaped_omegas = reshape(out.omegas', [nTrials * nBlocks, 4]);
            % else
            % end
            reshaped_omegas = reshape(out.omegas', [nTrials * nBlocks, 1]);
            if length(omegas) < i
                omegas{i} = reshaped_omegas; % Initialize the cell if it doesn't exist
            else
                omegas{i} = [omegas{i} reshaped_omegas]; % Append to the existing cell
            end
        end

        if isfield(out, 'omegasStim') && subj.selected_schedule_idx == 1
            reshaped_omegas = reshape(out.omegasStim, [nTrials * nBlocks, 4]);
            schedule_counter = schedule_counter + 1;


            omegasStim = omegasStim + reshaped_omegas;
        end
        
        if isfield(out, 'qs') && subj.selected_schedule_idx == 1
            reshaped_qs = reshape(out.qs', [nTrials * nBlocks, 1]);
            reshaped_svs = reshape(out.svs', [nTrials * nBlocks, 1]);

            if length(qs) < i
                qs{i} = reshaped_qs;
                svs{i} = reshaped_svs;
            else
                qs{i} = [qs{i} reshaped_qs];
                svs{i} = [svs{i} reshaped_svs];
            end
        end

      % Directly store results in the preallocated arrays
        allHCcell{iRun} = out.HCcell;
        allLCcell{iRun} = out.LCcell;
        allYCcell{iRun} = out.YCcell;

        allStayAnalysis(iRun, :, :) = calcSummaryStayAnalysis(out);
    end
    
    if isfield(out, 'omegasStim')
        allOmegaStim{k} = omegasStim./schedule_counter;
        omegasStim = zeros(nTrials * nBlocks, 4);
    end
    schedule_counter = 0;

    % Concatenate all results after the loop
    HCcell = vertcat(allHCcell{:});
    LCcell = vertcat(allLCcell{:});
    YCcell = vertcat(allYCcell{:});
    
    % Initialize variables to store sum and sum of squares
    sumStayAnalysis = zeros(size(allStayAnalysis, 2), size(allStayAnalysis, 3));
    sumSqStayAnalysis = zeros(size(allStayAnalysis, 2), size(allStayAnalysis, 3));
    
    % Loop over runs to calculate the sum and sum of squares
    for iRun = 1:nRuns
        sumStayAnalysis = sumStayAnalysis + squeeze(allStayAnalysis(iRun, :, :));
        sumSqStayAnalysis = sumSqStayAnalysis + squeeze(allStayAnalysis(iRun, :, :)).^2;
    end
    
    % Calculate the mean
    avgStayAnalysis = sumStayAnalysis / nRuns;
    
    % Calculate the standard error (assuming normal distribution)
    stderrStayAnalysis = sqrt((sumSqStayAnalysis - (sumStayAnalysis.^2 / nRuns)) / (nRuns - 1));
    stderrStayAnalysis = stderrStayAnalysis / sqrt(nRuns);  % Standard error
    

    % Plot the learning curves for each model
    figure(fig1);
    subplot(3, ceil(length(selMods)/3), i);
    plotLearningCurves(HCcell, sprintf("M%02d", idx), fig1);
    
    figure(fig2);
    subplot(3, ceil(length(selMods)/3), i);
    plotLearningCurves(LCcell, sprintf("M%02d", idx), fig2);
    
    figure(fig3)
    subplot(3, ceil(length(selMods)/3), i);
    plotLearningCurves(YCcell, sprintf("M%02d", idx), fig3);

    figure(fig4)
    subplot(3, ceil(length(selMods)/3), i);
    plotSummaryStayAnalysis(avgStayAnalysis, stderrStayAnalysis, sprintf("M%02d", idx), fig4);
end

% Initialize an empty cell array for average omegas
averageOmegas = {};
averageQs = {};
averageSvs = {};

% Counter to track the index in averageOmegas
index = 1;

% Loop through omegas to calculate the average for non-empty cells
for i = 1:length(omegas)
    if ~isempty(omegas{i})
        averageOmegas{index} = mean(omegas{i}, 2); % Average across columns
        index = index + 1;
    end
end

for i = 1:length(qs)
    if ~isempty(qs{i})
        averageQs{i} = mean(qs{i}, 2);
        averageSvs{i} = mean(svs{i}, 2);
    end
end
%% 02c) Plotting

% Plot the average Omegas for the first schedule
if ~isempty(omegas)
    plotOmegas(averageOmegas, cPs, modelsWithOmegas)
end

% Plot the average Omegas for the first schedule in stimulus indepedent
if any(omegasStim)
    plotOmegasStim(allOmegaStim, cPs, modelsWithOmegasStim)
end

if ~isempty(averageQs)
    plotQsandSvs(averageQs, averageSvs, cPs, modelsWithQs)
end
% % Plot the average Omegas
% figure(fig10);
% plotOmegas([averageOmegas6, averageOmegas7, averageOmegas8, averageOmegas9, averageOmegas10, averageOmegas11, averageOmegas12, averageOmegas13, averageOmegas14, averageOmegas15, averageOmegas16, averageOmegas17, averageOmegas18, averageOmegas19, averageOmegas20, averageOmegas21], cPs, fig10);

% Plot the average Omegas for the first schedule
% figure;
% hold on;
% colors = lines(nModelsWithOmega); % Get distinct colors for plotting
% for i = 1:length(averageOmegas)
%     if ~isempty(averageOmegas{i})
%         plot(averageOmegas{i}, 'Color', colors(i,:), 'DisplayName', sprintf('Model %d', modelsWithOmega(i)));
%     end
% end
% hold off;
% legend show;
% title('Average Omegas for Selected Models on the First Schedule');
% xlabel('Trial');
% ylabel('Omega value');
%% 3 Model inventor
% Quick turns between inventing, fitting and plotting. Meant to help with
% inventing new models quicker. 

% % Get separate datasets
% hc_d = cell(1, nSub);
% lc_d = cell(1, nSub);
% nBlocks = size(data{1}.controllability, 1);
% 
% for iSub = 1:nSub
%     d = data{iSub};
%     % Identify high control (hc) blocks
%     hc_idx = find(d.controllability(:, 1) == 1);
% 
%     % Identify low control (lc) blocks where controllability is 0 and isYoked is 0
%     lc_idx = find(d.controllability(:, 1) == 0);
% 
%     % Extract the data for the high control blocks and store it
%     hc_d{iSub} = structfun(@(x) x(hc_idx, :), d, 'UniformOutput', false);
%     lc_d{iSub} = structfun(@(x) x(lc_idx, :), d, 'UniformOutput', false);
% end

%% Select model and simulate it
% selMod = 22;
% modelHandle = str2func(sprintf("ChemControl_mod%d", selMod));
% modelSimHandle = str2func(sprintf("ChemControl_mod%d_modSim", selMod));
% 
% % Run settings
% nRuns = 1;
% nTrials = 80;
% nBlocks = 16;
% nStates = 4;
% nSchedules = 11;
% 
% % Simulate model with given parameters
% fname_mod = fullfile(dirs.lap, sprintf('lap_mod%02d_allData.mat', selMod));
% 
% % Check if the file exists
% if isfile(fname_mod) 
%     fname = load(fname_mod);
% else
%     % Fit model using LAP cbm if file does not exist
%     prior = struct('mean', [0 2 0 0 2 0 0 0], 'variance', [3 5 10 3 5 3 3 3]); % prior_model_dynamicOmega2
%     control_types = {'data', 'hc_d', 'lc_d'};
%     for ctype = control_types
%         ctype = ctype{1};
% 
%         fprintf('Processing %s data\n', ctype);
%         d = eval(ctype);
% 
%         fprintf('Fitting Laplace Approximation for %s data\n', ctype);
%         switch ctype
%             case 'data'
%                 fname_mod = fullfile(dirs.lap, sprintf('lap_mod%02d_allData.mat', selMod));
%             case 'hc_d'
%                 fname_mod = fullfile(dirs.lap, sprintf('lap_mod%02d_hc.mat', selMod));
%             case 'lc_d'
%                 fname_mod = fullfile(dirs.lap, sprintf('lap_mod%02d_lc.mat', selMod));
%         end
% 
%         cbm_lap(d, modelHandle, prior, fname_mod);
%     end
% end
% fname = load(fullfile(dirs.lap, sprintf('lap_mod%02d_allData.mat', selMod)));
% 
% cbm = fname.cbm;
% subParam  = cbm.output.parameters;
% parameters = mean(subParam, 1);
% 
% % Preallocate cell arrays for all outputs
% allHCcell = cell(nRuns, 1);
% allLCcell = cell(nRuns, 1);
% allYCcell = cell(nRuns, 1);
% 
% allOmegas = {}; % Store omega values for each run
% 
% for iRun = 1:nRuns
%     subj = sim_subj(nBlocks, nTrials);
% 
%     out = modelSimHandle(parameters, subj);
% 
%     % Directly store results in the preallocated arrays
%     allHCcell{iRun} = out.HCcell;
%     allLCcell{iRun} = out.LCcell;
%     allYCcell{iRun} = out.YCcell;
% 
%     % Store reshaped omega values
%     if subj.selected_schedule_idx == 1
%         allOmegas{end+1} = reshape(out.omegas', [nTrials*nBlocks, 1]);
%     end
% end
% % Concatenate all results after the loop
% HCcell = vertcat(allHCcell{:});
% LCcell = vertcat(allLCcell{:});
% YCcell = vertcat(allYCcell{:});
% 
% meanOmega = mean(allOmegas, 2);
% 
% figX = figure;
% sgtitle(sprintf("M%02d: Learning curves in different control types", selMod))
% 
% controls = ["HC", "LC", "YC"];
% for i = 1:length(controls)
%     subplot(1, 3, i)
%     plotLearningCurves(eval(sprintf("%scell", controls(i))), sprintf("%s control type", controls(i)), figX);
% end
% 
% figX = figure;
% sgtitle(sprintf("M%02d: Average Omega over time with controllability probability", selMod));
% hold on;
% plot(meanOmega, 'LineWidth', 2);
% plot(1:length(cPs(iSchedule, :)), cPs(iSchedule, :), 'LineWidth', 2, 'LineStyle', '--', 'Color', '#AEAEAE', 'DisplayName', 'Outcome Probability')
% grid on;
% hold off;

