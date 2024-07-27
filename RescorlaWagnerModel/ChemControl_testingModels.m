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
for iMod = 1:nMod
    % load
    fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d_%s.mat', iMod, dataType));
    fprintf("Load model %d fit with HBI\n", iMod);
    load(fname_hbi)
    fprintf("Extract model %d fit with HBI\n", iMod);
    groupParams{iMod} = cbm.output.group_mean{:};
end

%% 02) SIMULATE
%% Schedules

nRuns = 100;
nTrials = 80;
nBlocks = 16;
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

%% 02a) Run settings

fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig1)
sgtitle(sprintf("Learning curves for different models in high control trials for %i runs, %i blocks, %i trials", nRuns, nBlocks, nTrials))

fig2 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig2)
sgtitle(sprintf("Learning curves for different models in low control trials for %i runs, %i blocks, %i trials", nRuns, nBlocks, nTrials))

fig3 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig3)
sgtitle(sprintf("Learning curves for different models in yoked low control trials for %i runs, %i blocks, %i trials", nRuns, nBlocks, nTrials))

% fig4 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig4)
% sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins in High Control")
% 
% fig5 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig5);
% sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins in Low Non-Yoked Control")
% 
% fig6 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig6);
% sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins in Low Yoked Control")
% 
% fig7 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig7)
% sgtitle(sprintf("Prediction errors for different models in high control trials for %i runs", nRuns))
% 
% fig8 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig8)
% sgtitle(sprintf("Prediction errors for different models in low control trials for %i runs", nRuns))
% 
% fig9 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig9)
% sgtitle(sprintf("Prediction errors for different models in yoked low control trials for %i runs", nRuns))

fig10 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig10);

% fig11 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig11)
% sgtitle(sprintf("Average reward rate in high control trials for %i runs", nRuns))
% 
% fig12 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig12)
% sgtitle(sprintf("Average reward rate in low control trials for %i runs", nRuns))
% 
% fig13 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
% figure(fig13)
% sgtitle(sprintf("Average reward rate in yoked low control trials for %i runs", nRuns))

% Preallocate cell arrays for all outputs
allHCcell = cell(nRuns, 1);
allLCcell = cell(nRuns, 1);
allYCcell = cell(nRuns, 1);

omegas6 = cell(nSchedules, 1);
omegas7 = cell(nSchedules, 1);
omegas8 = cell(nSchedules, 1);
omegas9 = cell(nSchedules, 1);
omegas10 = cell(nSchedules, 1);
omegas11 = cell(nSchedules, 1);
omegas12 = cell(nSchedules, 1);
omegas13 = cell(nSchedules, 1);
omegas14 = cell(nSchedules, 1);
omegas15 = cell(nSchedules, 1);
omegas16 = cell(nSchedules, 1);
omegas17 = cell(nSchedules, 1);
omegas18 = cell(nSchedules, 1);
omegas19 = cell(nSchedules, 1);
omegas20 = cell(nSchedules, 1);
omegas21 = cell(nSchedules, 1);



selMods = 1:nMod;

% 02b) Run simulation and plot
for iMod = selMods

    parameters = groupParams{iMod};
    for iRun = 1:nRuns
        subj = sim_subj(nBlocks, nTrials);
        out = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", iMod));
        % if iMod >= 7 && iMod <= 21
        %     % Store omegas for models 6 to 21
        %     schedule_idx = subj.selected_schedule_idx;
        %     reshaped_omegas = reshape(out.omegas', [nTrials*nBlocks, 1]);
        % 
        %     if iMod == 6
        %         if isempty(omegas6{schedule_idx})
        %             omegas6{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas6{schedule_idx} = [omegas6{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 7
        %         if isempty(omegas7{schedule_idx})
        %             omegas7{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas7{schedule_idx} = [omegas7{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 8
        %         if isempty(omegas8{schedule_idx})
        %             omegas8{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas8{schedule_idx} = [omegas8{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 9
        %         if isempty(omegas9{schedule_idx})
        %             omegas9{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas9{schedule_idx} = [omegas9{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 10
        %         if isempty(omegas10{schedule_idx})
        %             omegas10{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas10{schedule_idx} = [omegas10{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 11
        %         if isempty(omegas11{schedule_idx})
        %             omegas11{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas11{schedule_idx} = [omegas11{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 12
        %         if isempty(omegas12{schedule_idx})
        %             omegas12{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas12{schedule_idx} = [omegas12{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 13
        %         if isempty(omegas13{schedule_idx})
        %             omegas13{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas13{schedule_idx} = [omegas13{schedule_idx}, reshaped_omegas];
        %         end            
        %     elseif iMod == 14
        %         if isempty(omegas14{schedule_idx})
        %             omegas14{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas14{schedule_idx} = [omegas14{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 15
        %         if isempty(omegas15{schedule_idx})
        %             omegas15{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas15{schedule_idx} = [omegas15{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 16
        %         if isempty(omegas16{schedule_idx})
        %             omegas16{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas16{schedule_idx} = [omegas16{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 17
        %         if isempty(omegas17{schedule_idx})
        %             omegas17{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas17{schedule_idx} = [omegas17{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 18
        %         if isempty(omegas18{schedule_idx})
        %             omegas18{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas18{schedule_idx} = [omegas18{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 19
        %         if isempty(omegas19{schedule_idx})
        %             omegas19{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas19{schedule_idx} = [omegas19{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 20
        %         if isempty(omegas20{schedule_idx})
        %             omegas20{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas20{schedule_idx} = [omegas20{schedule_idx}, reshaped_omegas];
        %         end
        %     elseif iMod == 21
        %         if isempty(omegas21{schedule_idx})
        %             omegas21{schedule_idx} = reshaped_omegas;
        %         else
        %             omegas21{schedule_idx} = [omegas21{schedule_idx}, reshaped_omegas];
        %         end
        %     end
        % end

      % Directly store results in the preallocated arrays
        allHCcell{iRun} = out.HCcell;
        allLCcell{iRun} = out.LCcell;
        allYCcell{iRun} = out.YCcell;
    
    end

    % Concatenate all results after the loop
    HCcell = vertcat(allHCcell{:});
    LCcell = vertcat(allLCcell{:});
    YCcell = vertcat(allYCcell{:});
  
    figure(fig1);
    subplot(3, ceil(nMod/3), iMod);
    plotLearningCurves(HCcell, sprintf("M%02d", iMod), fig1);
    
    figure(fig2);
    subplot(3, ceil(nMod/3), iMod);
    plotLearningCurves(LCcell, sprintf("M%02d", iMod), fig2);
    
    figure(fig3)
    subplot(3, ceil(nMod/3), iMod);
    plotLearningCurves(YCcell, sprintf("M%02d", iMod), fig3);
end

% Initialize cell arrays for average Omegas
averageOmegas6 = cell(nSchedules, 1);
averageOmegas7 = cell(nSchedules, 1);
averageOmegas8 = cell(nSchedules, 1);
averageOmegas9 = cell(nSchedules, 1);
averageOmegas10 = cell(nSchedules, 1);
averageOmegas11 = cell(nSchedules, 1);
averageOmegas12 = cell(nSchedules, 1);
averageOmegas13 = cell(nSchedules, 1);
averageOmegas14 = cell(nSchedules, 1);
averageOmegas15 = cell(nSchedules, 1);
averageOmegas16 = cell(nSchedules, 1);
averageOmegas17 = cell(nSchedules, 1);
averageOmegas18 = cell(nSchedules, 1);
averageOmegas19 = cell(nSchedules, 1);
averageOmegas20 = cell(nSchedules, 1);
averageOmegas21 = cell(nSchedules, 1);

for schedule_idx = 1:nSchedules
    if ~isempty(omegas6{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas6{schedule_idx} = mean(omegas6{schedule_idx}, 2);
    end
    if ~isempty(omegas7{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas7{schedule_idx} = mean(omegas7{schedule_idx}, 2);
    end
    if ~isempty(omegas8{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas8{schedule_idx} = mean(omegas8{schedule_idx}, 2);
    end
    if ~isempty(omegas9{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas9{schedule_idx} = mean(omegas9{schedule_idx}, 2);
    end
    if ~isempty(omegas10{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas10{schedule_idx} = mean(omegas10{schedule_idx}, 2);
    end
    if ~isempty(omegas11{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas11{schedule_idx} = mean(omegas11{schedule_idx}, 2);
    end
    if ~isempty(omegas12{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas12{schedule_idx} = mean(omegas12{schedule_idx}, 2);
    end
    if ~isempty(omegas13{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas13{schedule_idx} = mean(omegas13{schedule_idx}, 2);
    end
    if ~isempty(omegas14{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas14{schedule_idx} = mean(omegas14{schedule_idx}, 2);
    end
    if ~isempty(omegas15{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas15{schedule_idx} = mean(omegas15{schedule_idx}, 2);
    end
    if ~isempty(omegas16{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas16{schedule_idx} = mean(omegas16{schedule_idx}, 2);
    end
    if ~isempty(omegas17{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas17{schedule_idx} = mean(omegas17{schedule_idx}, 2);
    end
    if ~isempty(omegas18{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas18{schedule_idx} = mean(omegas18{schedule_idx}, 2);
    end
    if ~isempty(omegas19{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas19{schedule_idx} = mean(omegas19{schedule_idx}, 2);
    end
    if ~isempty(omegas20{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas20{schedule_idx} = mean(omegas20{schedule_idx}, 2);
    end
    if ~isempty(omegas21{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas21{schedule_idx} = mean(omegas21{schedule_idx}, 2);
    end
end

% Plot the average Omegas
figure(fig10);
plotOmegas([averageOmegas6, averageOmegas7, averageOmegas8, averageOmegas9, averageOmegas10, averageOmegas11, averageOmegas12, averageOmegas13, averageOmegas14, averageOmegas15, averageOmegas16, averageOmegas17, averageOmegas18, averageOmegas19, averageOmegas20, averageOmegas21], cPs, fig10);

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

