% ChemControl_analyzeData.m
% 
% Execute this script to run analyses for the preparedData.
% ChemControl_cbm_prepareData.m has to be run before.
%
% INPUTS:
% none.
%
% OUTPUTS:
% No outputs, just plots
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/

clear all; close all; clc

%% Directories

dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

%% Input file

inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
data = load(inputFile).data;
nSub = size(data, 2);
fprintf('Loaded inputFile with %d subjects\n', size(data, 2));

%% Settings
nSub = size(data,2);
nBlocks = size(data{1}.actions, 1);
nTrials = size(data{1}.actions, 2);

%% PLOTS
% Learning curves
fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig1)
sgtitle(sprintf("Learning Curves for Subjects in High Control Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

fig2 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig2)
sgtitle(sprintf("Learning Curves for Subjects in Low, Non-yoked Control Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

fig3 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig3)
sgtitle(sprintf("Learning Curves for Subjects in Low, Yoked Control Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

fig4 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig4)
sgtitle(sprintf("Learning Curves for Subjects in All Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

%% Extracting subsets of data
hc_d = cell(1, nSub);
lc_d = cell(1, nSub);
yc_d = cell(1, nSub);
nBlocks = size(data{1}.controllability, 1);
for iSub = 1:nSub
    d = data{iSub};  % Get the struct for the current subject
    
    % Identify high control (hc) blocks
    hc_idx = find(d.controllability(:, 1));
    
    % Identify low control (lc) blocks where controllability is 0 and isYoked is 0
    lc_idx = find(~d.controllability(:, 1) & ~d.isYoked(:, 1));
    
    % Identify yoked control (yc) blocks where isYoked is 1
    yc_idx = find(d.isYoked(:, 1));
    
    % Extract the data for the high control blocks and store it
    hc_d{iSub} = structfun(@(x) x(hc_idx, :), d, 'UniformOutput', false);
    lc_d{iSub} = structfun(@(x) x(lc_idx, :), d, 'UniformOutput', false);
    yc_d{iSub} = structfun(@(x) x(yc_idx, :), d, 'UniformOutput', false);
end

%% Plot learning curves for each data type

plotParticipantCurves(hc_d, fig1);
plotParticipantCurves(lc_d, fig2);
plotParticipantCurves(yc_d, fig3);
plotParticipantCurves(data, fig4);

%% Block transitions with conficence bounds
% Initialize data structures for each transition
transitions = {'hchc', 'hclc', 'hcyc', 'lclc', 'lchc', 'lcyc', 'ycyc', 'yclc', 'ychc'};
transitionData = struct();

for i = 1:numel(transitions)
    transitionData.(transitions{i}) = {};
end

% Loop over subjects
for iSub = 1:nSub
    subj = data{iSub};
    schedule = subj.controllability(:, 1) + 2* subj.isYoked(:, 1);
    for iBlock = 1:(numel(schedule) - 1)
        % Convert transition to a string key
        transition = sprintf('%d%d', schedule(iBlock), schedule(iBlock + 1));

        switch transition
            case '11'
                key = 'hchc';
            case '10'
                key = 'hclc';
            case '12'
                key = 'hcyc';
            case '00'
                key = 'lclc';
            case '01'
                key = 'lchc';
            case '02'
                key = 'lcyc';
            case '22'
                key = 'ycyc';
            case '20'
                key = 'yclc';
            case '21'
                key = 'ychc';
            otherwise
                continue;
        end

        % Collect the data for the current transition
        transitionData.(key){end + 1} = structfun(@(x) x(iBlock + 1, :), subj, 'UniformOutput', false);
    end
end

fig5 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig5)
sgtitle(sprintf("Learning Curves for Specific Transition Types for Different Data Types for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

plotTransitions(transitionData, fig5);




