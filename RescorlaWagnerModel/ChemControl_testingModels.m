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
nRuns = 50;
nTrials = 40;
nBlocks = 32;
nStates = 4;
nSchedules = 11;
HCmeans = zeros(nTrials/4, 4, nRuns);
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
fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig1)
sgtitle("Learning curves for different models in high control trials")

fig2 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
fig3 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig3)
sgtitle("Weighted Probability of Shifting after a Loss vs. Number of Consecutive Wins")



% 02b) Run simulation and plot
for iMod = 1:nMod
    shiftMeans = zeros(nRuns, nStates, nTrials);
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
        HCoccurrences = zeros(nTrials/4, 4);
        for s = 1:4
            for occurrence = 1:nTrials/4
                HCprobs = zeros(nBlocks/2, 1);
                for b = 1:nBlocks/2
                    HCprobs(b) = HCcell{b, s}(occurrence);
                end
                HCoccurrences(occurrence, s) = mean(HCprobs);
            end
        end
        HCmeans(:, :, iRun) = HCoccurrences; 
        shiftMeans(iRun, :, :) = out.weightedProbShiftAfterLoss;

    end
    HCmeans = mean(HCmeans, 3);
    figure(fig1);
    subplot(2, ceil(nMod/2), iMod);
    plotLearningCurves(HCmeans, sprintf("M%02d", iMod), fig1);

    shiftMeans = squeeze(mean(shiftMeans, 1));
    figure(fig3);
    subplot(2, ceil(nMod/2), iMod);
    plotShiftLoss(shiftMeans, iMod, fig3);
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

figure(fig2);
plotOmegas(averageOmegas6, averageOmegas7, cPs, fig2);




