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
dirs.models         = fullfile(dirs.root, 'models');
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

% Input file
inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
data = load(inputFile).data;
fprintf('Loaded inputFile with %d subjects\n', size(data, 2));

nMod = length(dir(fullfile(dirs.models, "*.m")));

%% 01) LOAD AND SAVE:
%% 01a) Alternative A: Load model fitted with LAP:

parType = 'lap';
selMod = 1;
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
% extract
for iMod = 1:nMod
    % load
    fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d.mat', iMod));
    fprintf("Load model %d fit with HBI\n", iMod);
    load(fname_hbi)
    fprintf("Extract model %d fit with HBI\n", iMod);
    groupParams{iMod} = cbm.output.group_mean{:};
end

%% 01c) Alternative C: Set parameters manually:
...

%% 02) SIMULATE
%% 02a) Run settings

nRuns = 100;
nTrials = 80;
nBlocks = 32;
nStates = 4;
nSchedules = 11;

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
sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins in High Control")

fig5 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig5);
sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins in Low Non-Yoked Control")

fig6 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig6);
sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins in Low Yoked Control")

fig7 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig7)
sgtitle(sprintf("Prediction errors for different models in high control trials for %i runs", nRuns))

fig8 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig8)
sgtitle(sprintf("Prediction errors for different models in low control trials for %i runs", nRuns))

fig9 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig9)
sgtitle(sprintf("Prediction errors for different models in yoked low control trials for %i runs", nRuns))

fig10 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig10);

fig11 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig11)
sgtitle(sprintf("Average reward rate in high control trials for %i runs", nRuns))

fig12 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig12)
sgtitle(sprintf("Average reward rate in low control trials for %i runs", nRuns))

fig13 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig13)
sgtitle(sprintf("Average reward rate in yoked low control trials for %i runs", nRuns))

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

HCmeans = zeros(nTrials/4, 4, nRuns);
LCmeans = zeros(nTrials/4, 4, nRuns);
YCmeans = zeros(nTrials/4, 4, nRuns);

HCpes_mean = zeros(nTrials/4, 4, nRuns);
LCpes_mean = zeros(nTrials/4, 4, nRuns);
YCpes_mean = zeros(nTrials/4, 4, nRuns);

HCarr_mean = zeros(nTrials/4, 4, nRuns);
LCarr_mean = zeros(nTrials/4, 4, nRuns);
YCarr_mean = zeros(nTrials/4, 4, nRuns);

omegas6 = cell(nSchedules, 1);
omegas7 = cell(nSchedules, 1);
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

% 02b) Run simulation and plot
for iMod = 1:nMod
    shiftMeans_HC = zeros(nRuns, nStates, nTrials);
    shiftMeans_LC = zeros(nRuns, nStates, nTrials);
    shiftMeans_YC = zeros(nRuns, nStates, nTrials);
    parameters = groupParams{iMod};
    for iRun = 1:nRuns
        subj = sim_subj(nBlocks, nTrials);
        out = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", iMod));
        if iMod == 6 || iMod == 7
            % Store omegas for models 6 and 7
            schedule_idx = subj.selected_schedule_idx;
            reshaped_omegas = reshape(out.omegas', [nTrials*nBlocks, 1]);
            
            if iMod == 6
                if isempty(omegas6{schedule_idx})
                    omegas6{schedule_idx} = reshaped_omegas;
                else
                    omegas6{schedule_idx} = [omegas6{schedule_idx}, reshaped_omegas];
                end
            elseif iMod == 7
                if isempty(omegas7{schedule_idx})
                    omegas7{schedule_idx} = reshaped_omegas;
                else
                    omegas7{schedule_idx} = [omegas7{schedule_idx}, reshaped_omegas];
                end            
            end
        end
        HCcell = out.HCcell;
        LCcell = out.LCcell;
        YCcell = out.YCcell;
        
        HCpes = out.HCpe;
        LCpes = out.LCpe;
        YCpes = out.YCpe;

        HCarr = out.HCarr;
        LCarr = out.LCarr;
        YCarr = out.YCarr;

        HCoccurrences = zeros(nTrials/4, 4);
        LCoccurrences = zeros(nTrials/4, 4);
        YCoccurrences = zeros(nTrials/4, 4);

        HCpe_occurrence = zeros(nTrials/4, 4);
        LCpe_occurrence = zeros(nTrials/4, 4);
        YCpe_occurrence = zeros(nTrials/4, 4);

        HCarr_occurrence = zeros(nTrials/4, 4);
        LCarr_occurrence = zeros(nTrials/4, 4);
        YCarr_occurrence = zeros(nTrials/4, 4);

        for s = 1:4
            for occurrence = 1:nTrials/4
                HCprobs = zeros(nBlocks/2, 1);
                LCprobs = zeros(nBlocks/4, 1);
                YCprobs = zeros(nBlocks/4, 1);

                HCpes_block = zeros(nBlocks/2, 1);
                LCpes_block = zeros(nBlocks/4, 1);
                YCpes_block = zeros(nBlocks/4, 1);

                HCarr_block = zeros(nBlocks/2, 1);
                LCarr_block = zeros(nBlocks/4, 1);
                YCarr_block = zeros(nBlocks/4, 1);

                for b = 1:nBlocks/2
                    HCprobs(b) = HCcell{b, s}(occurrence);
                    HCpes_block(b) = HCpes{b, s}(occurrence);
                    HCarr_block(b) = HCarr{b, s}(occurrence);
                    
                    if b <= nBlocks/4
                        LCprobs(b) = LCcell{b, s}(occurrence);
                        YCprobs(b) = YCcell{b, s}(occurrence);

                        LCpes_block(b) = LCpes{b, s}(occurrence);
                        YCpes_block(b) = YCpes{b, s}(occurrence);

                        LCarr_block(b) = LCarr{b, s}(occurrence);
                        YCarr_block(b) = YCarr{b, s}(occurrence);
                    end
                end
                HCoccurrences(occurrence, s) = mean(HCprobs);
                LCoccurrences(occurrence, s) = mean(LCprobs);
                YCoccurrences(occurrence, s) = mean(YCprobs);

                HCpe_occurrence(occurrence, s) = mean(HCpes_block);
                LCpe_occurrence(occurrence, s) = mean(LCpes_block);
                YCpe_occurrence(occurrence, s) = mean(YCpes_block);

                HCarr_occurrence(occurrence, s) = mean(HCarr_block);
                LCarr_occurrence(occurrence, s) = mean(LCarr_block);
                YCarr_occurrence(occurrence, s) = mean(YCarr_block);
           end
        end
        HCmeans(:, :, iRun) = HCoccurrences; 
        LCmeans(:, :, iRun) = LCoccurrences;
        YCmeans(:, :, iRun) = YCoccurrences;

        HCpes_mean(:, :, iRun) = HCpe_occurrence;
        LCpes_mean(:, :, iRun) = LCpe_occurrence;
        YCpes_mean(:, :, iRun) = YCpe_occurrence;

        HCarr_mean(:, :, iRun) = HCarr_occurrence; 
        LCarr_mean(:, :, iRun) = LCarr_occurrence;
        YCarr_mean(:, :, iRun) = YCarr_occurrence;

        % Store shift means for each control type
        shiftMeans_HC(iRun, :, :) = out.weightedProbShiftAfterLoss_HC;
        shiftMeans_LC(iRun, :, :) = out.weightedProbShiftAfterLoss_LC;
        shiftMeans_YC(iRun, :, :) = out.weightedProbShiftAfterLoss_YC;

    end
    HCmeans = mean(HCmeans, 3);
    LCmeans = mean(LCmeans, 3);
    YCmeans = mean(YCmeans, 3);

    HCpes_mean = mean(HCpes_mean, 3);
    LCpes_mean = mean(LCpes_mean, 3);
    YCpes_mean = mean(YCpes_mean, 3);

    HCarr_mean = mean(HCarr_mean, 3);
    LCarr_mean = mean(LCarr_mean, 3);
    YCarr_mean = mean(YCarr_mean, 3);

    figure(fig1);
    subplot(2, ceil(nMod/2), iMod);
    plotLearningCurves(HCmeans, sprintf("M%02d", iMod), fig1);
    
    figure(fig2);
    subplot(2, ceil(nMod/2), iMod);
    plotLearningCurves(LCmeans, sprintf("M%02d", iMod), fig2);
    
    figure(fig3)
    subplot(2, ceil(nMod/2), iMod);
    plotLearningCurves(YCmeans, sprintf("M%02d", iMod), fig3);

    shiftMeans_HC = squeeze(mean(shiftMeans_HC, 1));
    shiftMeans_LC = squeeze(mean(shiftMeans_LC, 1));
    shiftMeans_YC = squeeze(mean(shiftMeans_YC, 1));

    figure(fig4);
    subplot(2, ceil(nMod/2), iMod);
    plotShiftLoss(shiftMeans_HC, iMod, fig4);

    figure(fig5);
    subplot(2, ceil(nMod/2), iMod);
    plotShiftLoss(shiftMeans_LC, iMod, fig5);

    figure(fig6);
    subplot(2, ceil(nMod/2), iMod);
    plotShiftLoss(shiftMeans_YC, iMod, fig6);

    figure(fig7);
    subplot(2, ceil(nMod/2), iMod);
    plotPredictionErrors(HCpes_mean, iMod, fig7);

    figure(fig8);
    subplot(2, ceil(nMod/2), iMod);
    plotPredictionErrors(LCpes_mean, iMod, fig8);

    figure(fig9);
    subplot(2, ceil(nMod/2), iMod);
    plotPredictionErrors(YCpes_mean, iMod, fig9);

    figure(fig11);
    subplot(2, ceil(nMod/2), iMod);
    plotAverageRewardRate(HCarr_mean, sprintf("M%02d", iMod), fig11);
    
    figure(fig12);
    subplot(2, ceil(nMod/2), iMod);
    plotAverageRewardRate(LCarr_mean, sprintf("M%02d", iMod), fig12);
    
    figure(fig13)
    subplot(2, ceil(nMod/2), iMod);
    plotAverageRewardRate(YCarr_mean, sprintf("M%02d", iMod), fig13);
end

% Plot average Omega
averageOmegas6 = cell(nSchedules, 1);
averageOmegas7 = cell(nSchedules, 1);


for schedule_idx = 1:nSchedules
    if ~isempty(omegas6{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas6{schedule_idx} = mean(omegas6{schedule_idx}, 2);
    end
    if ~isempty(omegas7{schedule_idx})
        % Average across columns for schedule_idx
        averageOmegas7{schedule_idx} = mean(omegas7{schedule_idx}, 2);
    end
end

figure(fig10);
plotOmegas(averageOmegas6, averageOmegas7, cPs, fig10);




