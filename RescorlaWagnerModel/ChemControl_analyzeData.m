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
dirs.hbi = fullfile(dirs.target, 'HBI_Results');

%% Input file

inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
data = load(inputFile).data;
nSub = size(data, 2);
fprintf('Loaded inputFile with %d subjects\n', size(data, 2));

inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData_invalidSubs.mat');
invalidData = load(inputFile).data;
invalidNSub = size(invalidData, 2);
fprintf('Loaded invalid inputFile with %d subjects\n', size(invalidData, 2));

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

fig8 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig8)
sgtitle(sprintf("Learning Curves for invalid subjects in High Control Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

fig9 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig9)
sgtitle(sprintf("Learning Curves for invalid subjects in Low, Non-yoked Control Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

fig10 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig10)
sgtitle(sprintf("Learning Curves for invalid subjects in Low, Yoked Control Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

fig11 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig11)
sgtitle(sprintf("Learning Curves for invalid subjects in All Trials for %i Subs, %i Blocks, %i Trials", nSub, nBlocks, nTrials))

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

%% Extracting subsets of data for invalid subjects
hc_invalid_d = cell(1, invalidNSub);
lc_invalid_d = cell(1, invalidNSub);
yc_invalid_d = cell(1, invalidNSub);
nBlocksInvalid = size(invalidData{1}.controllability, 1);

for iSub = 1:invalidNSub
    d_invalid = invalidData{iSub};  % Get the struct for the current invalid subject
    
    % Identify high control (hc) blocks
    hc_idx_invalid = find(d_invalid.controllability(:, 1));
    
    % Identify low control (lc) blocks where controllability is 0 and isYoked is 0
    lc_idx_invalid = find(~d_invalid.controllability(:, 1) & ~d_invalid.isYoked(:, 1));
    
    % Identify yoked control (yc) blocks where isYoked is 1
    yc_idx_invalid = find(d_invalid.isYoked(:, 1));
    
    % Extract the data for the high control blocks and store it
    hc_invalid_d{iSub} = structfun(@(x) x(hc_idx_invalid, :), d_invalid, 'UniformOutput', false);
    lc_invalid_d{iSub} = structfun(@(x) x(lc_idx_invalid, :), d_invalid, 'UniformOutput', false);
    yc_invalid_d{iSub} = structfun(@(x) x(yc_idx_invalid, :), d_invalid, 'UniformOutput', false);
end

%% Plot learning curves for each data type for invalid subjects

plotParticipantCurves(hc_invalid_d, fig8);
plotParticipantCurves(lc_invalid_d, fig9);
plotParticipantCurves(yc_invalid_d, fig10);
plotParticipantCurves(invalidData, fig11);

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

%% Average reward rate with running windows
windowSize = 7;

fig6 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig6)
sgtitle(sprintf("Average Reward Rate for %i Subs, with window size of %i", nSub, windowSize))

ctypes = {"all", "hc", "lc", "yc"};
for i = 1:numel(ctypes)
    ctype = ctypes{i};

    if strcmp(ctype, "all")
        currentData = data;
    elseif strcmp(ctype, "hc")
        currentData = hc_d;
    elseif strcmp(ctype, "lc")
        currentData = lc_d;
    elseif strcmp(ctype, "yc")
        currentData = yc_d;
    end
    nBlocks = size(currentData{1}.outcomes, 1);
    nTrials = size(currentData{1}.outcomes, 2);

    allAverageRewardRates = zeros(nBlocks, nTrials - windowSize, nSub);

    for iSub = 1:nSub
        subj = currentData{iSub};
        outcomes = subj.outcomes;
        outcomes = transformOutcomes(outcomes, subj.stimuli);
        for iBlock = 1:nBlocks
            for iTrial = windowSize + 1:nTrials
                window = outcomes(iBlock, iTrial-windowSize:iTrial);
                allAverageRewardRates(iBlock, iTrial - windowSize, iSub) = mean(window);
            end
        end
    end
    
    meanAverageRewardRate = mean(allAverageRewardRates, 3);
    meanAverageRewardRate = mean(meanAverageRewardRate, 1);  % Mean across blocks

    % Plotting for each control type
    subplot(2, 2, i)
    plot(meanAverageRewardRate, 'LineWidth', 2);
    title(['Average Reward Rate - ', ctype]);
    xlabel('Trial (adjusted for window size)');
    ylabel('Average Reward Rate');
    xlim([1 length(meanAverageRewardRate)])
    ylim([0.5 0.85]);
    grid on;
end

%% Average reward rate across transitions with sliding window

transitionArr = struct();

% Initialize each transition type in the structure
for i = 1:numel(transitions)
    transitionArr.(transitions{i}).stimuli = [];
    transitionArr.(transitions{i}).outcomes = [];
end

% Loop over subjects
for iSub = 1:nSub
    subj = data{iSub};
    schedule = subj.controllability(:, 1) + 2 * subj.isYoked(:, 1);
    nBlocks = numel(schedule);

    for iBlock = 1:(nBlocks - 1)
        % Convert transition to a string key
        transition = sprintf('%d%d', schedule(iBlock), schedule(iBlock + 1));
        key = '';

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

        if isempty(key)
            continue;
        end

        % Define indices for window around the transition
        endIdx = size(subj.stimuli, 2);  % Number of trials in each block
        preStart = max(1, endIdx - windowSize*2 + 1);
        postEnd = min(endIdx, windowSize*2);

        % Collect data from the end of the current block and the beginning of the next block
        if iBlock < nBlocks
            stimuliWindow = [subj.stimuli(iBlock, preStart:end), subj.stimuli(iBlock + 1, 1:postEnd)];
            outcomesWindow = [subj.outcomes(iBlock, preStart:end), subj.outcomes(iBlock + 1, 1:postEnd)];
            transitionArr.(key).stimuli = [transitionArr.(key).stimuli; stimuliWindow];
            transitionArr.(key).outcomes = [transitionArr.(key).outcomes; outcomesWindow];
        end
    end
end

transitions = fieldnames(transitionArr);
fig7 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig7)
sgtitle(sprintf("Average reward rate transitions for %i Subs, with windowSize %i ", nSub, windowSize))

for idx = 1:numel(transitions)
    transition = transitions{idx};
    transitionOutcomes = transitionArr.(transition).outcomes;
    transitionStimuli = transitionArr.(transition).stimuli;
    
    transitionOutcomes = transformOutcomes(transitionOutcomes, transitionStimuli);

    nTransitions = size(transitionStimuli, 1);
    nTrialsTransition = size(transitionStimuli, 2);

    allTransitionArrs = zeros(nTransitions, nTrialsTransition - windowSize);
    for iTransition = 1:nTransitions
        for iTrial = windowSize + 1:nTrialsTransition
            window = transitionOutcomes(iTransition, iTrial-windowSize:iTrial);
            allTransitionArrs(iTransition, iTrial - windowSize) = mean(window);
        end
    end
    
    meanTransitionArr = mean(allTransitionArrs, 1);
    subplot(3, 3, idx)
    % Plotting the average reward rate for each transition type
    plot(meanTransitionArr, 'LineWidth',2);
    title(['Average Reward Rate for ', transition, ' Transition']);
    xlabel('Trials (adjusted for window size)');
    ylabel('Average Reward Rate');
    xlim([1 length(meanTransitionArr)])
    ylim([0.5 0.85]);
    grid on;
end

accBySub = zeros(nSub, 1);

for iSub = 1:nSub
    subj = data{iSub};

    actions = subj.actions;
    states = subj.stimuli;
    goStates = states <= 2; 
    noGoStates = states >= 2;

    actions = actions == 1;
    correctActions = goStates & actions | noGoStates & ~actions;
    acc = mean(correctActions(:, 30:40),'all');
    accBySub(iSub) = acc;
end

% Calculate mean and standard deviation
meanAcc = mean(accBySub);
stdAcc = std(accBySub);

% Identify subjects more than 2 standard deviations below the mean
threshold = meanAcc - 2 * stdAcc;
outliers = accBySub < threshold;

fig8 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
figure(fig8)
sgtitle("Accuracy by Subject")
plot(accBySub, 'LineWidth', 2);
hold on;
yline(0.5, "LineStyle", "-", "Color", "#808080", 'LineWidth', 2)

% Highlight outliers
plot(find(outliers), accBySub(outliers), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

ylabel("Accuracy")
xlabel("Subject")
grid on;

% Annotate outliers
for i = find(outliers)'
    text(i, accBySub(i), sprintf('%.2f', accBySub(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
%% Response time by state
allHcRts = {};
allLcRts = {};
allYcRts = {};

for iSub = 1:nSub
    subj = data{iSub};
    schedule = 2 * (~subj.controllability(:, 1) & subj.isYoked(:, 1)) + subj.controllability(:, 1);
    for iBlock = 1:nBlocks
        blockGos = subj.actions(iBlock, :) == 1;
        blockRts = subj.responseTime(iBlock, blockGos);
        switch schedule(iBlock)
            case 1
                allHcRts = vertcat(allHcRts, blockRts);
            case 0
                allLcRts = vertcat(allLcRts, blockRts);
            case 2
                allYcRts = vertcat(allYcRts, blockRts);
        end
    end
end



% Calculate mean response times for each condition
meanHcResponseTime = calculateMeanResponseTime(allHcRts);
meanLcResponseTime = calculateMeanResponseTime(allLcRts);
meanYcResponseTime = calculateMeanResponseTime(allYcRts);

% Display the results
fprintf('Mean HC Response Time: %.4f seconds\n', mean(meanHcResponseTime));
fprintf('Mean LC Response Time: %.4f seconds\n', mean(meanLcResponseTime));
fprintf('Mean YC Response Time: %.4f seconds\n', mean(meanYcResponseTime));

% Optional: Plot the results
figure;
hold on;
plot(meanHcResponseTime, 'LineWidth', 2, 'DisplayName', 'High Control');
plot(meanLcResponseTime, 'LineWidth', 2, 'DisplayName', 'Low Control');
plot(meanYcResponseTime, 'LineWidth', 2, 'DisplayName', 'Yoked Control');
xlabel('Trial (up to mean length)');
ylabel('Response Time (seconds)');
legend show;
title('Mean Response Times by Condition');
grid on;
hold off;




%% Initialize variables to track exploratory actions
allHcExploration = [];
allLcExploration = [];
allYcExploration = [];

for iSub = 1:nSub
    subj = data{iSub};
    schedule = 2 * (~subj.controllability(:, 1) & subj.isYoked(:, 1)) + subj.controllability(:, 1);
    
    for iBlock = 1:nBlocks
        blockActions = subj.actions(iBlock, :);
        blockOutcomes = subj.outcomes(iBlock, :);
        blockStimuli = subj.stimuli(iBlock, :);
        % Determine the optimal action based on previous outcomes
        optimalActions = zeros(size(blockActions));
        for t = 2:nTrials
            % Optimal action is the one that was most frequently rewarded in the past
            pastOutcomes = blockOutcomes(1:t-1);
            pastActions = blockActions(1:t-1);
            pastStimuli = blockStimuli(1:t-1);

            % Calculate the number of rewards for each action considering the stimuli state
            rewardsAction1 = sum((pastOutcomes == 1 & mod(pastStimuli, 2) == 1 & pastActions == 1) | (pastOutcomes == 0 & mod(pastStimuli, 2) == 0 & pastActions == 1));
            rewardsAction2 = sum((pastOutcomes == 1 & mod(pastStimuli, 2) == 1 & pastActions == 2) | (pastOutcomes == 0 & mod(pastStimuli, 2) == 0 & pastActions == 2));
            if rewardsAction1 >= rewardsAction2
                optimalActions(t) = 1;
            else
                optimalActions(t) = 2;
            end
        end
        
        % Calculate exploratory actions
        exploratoryActions = blockActions ~= optimalActions;
        exploratoryRate = mean(exploratoryActions);
        
        % Store exploratory rates based on control condition
        switch schedule(iBlock)
            case 1
                allHcExploration = [allHcExploration; exploratoryRate];
            case 0
                allLcExploration = [allLcExploration; exploratoryRate];
            case 2
                allYcExploration = [allYcExploration; exploratoryRate];
        end
    end
end

% Calculate the average exploration rate per trial
meanHcExplorationRate = mean(allHcExploration);
meanLcExplorationRate = mean(allLcExploration);
meanYcExplorationRate = mean(allYcExploration);

% Calculate the standard error of the mean (SEM) for each condition
semHcExplorationRate = std(allHcExploration) / sqrt(length(allHcExploration));
semLcExplorationRate = std(allLcExploration) / sqrt(length(allLcExploration));
semYcExplorationRate = std(allYcExploration) / sqrt(length(allYcExploration));

% Create a bar plot with error bars
figure;
hold on;
barData = [meanHcExplorationRate, meanLcExplorationRate, meanYcExplorationRate];
semData = [semHcExplorationRate, semLcExplorationRate, semYcExplorationRate];
b = bar(barData);
x = 1:length(barData);
errorbar(x, barData, semData, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTick', x, 'XTickLabel', {'High Control', 'Low Control', 'Yoked Control'});
ylabel('Mean Exploration Rate');
title('Exploration Rate by Control Condition');

% Add error bars
numGroups = size(barData, 2);
x = 1:numGroups;
errorbar(x, barData, semData, 'k', 'linestyle', 'none', 'LineWidth', 1);

% Customize the plot
grid on;
hold off;

%% Win Stay loose Shift Analysis

%% Initialize variables to track exploratory actions and outcomes for all subjects
weightedProbShiftAfterLoss_HC_perSubj = {};
weightedProbShiftAfterLoss_LC_perSubj = {};
weightedProbShiftAfterLoss_YC_perSubj = {};

for iSub = 1:nSub
    % Filtering data
    highControlStimuli = [];
    lowControlStimuli = [];
    yokedLowControlStimuli = [];
    
    highControlActions = [];
    lowControlActions = [];
    yokedLowControlActions = [];
    
    highControlOutcomes = [];
    lowControlOutcomes = [];
    yokedLowControlOutcomes = [];
    
    hc = 0;
    lc = 0;
    yc = 0; % Counter for YokedLowControl
    
    subj = data{iSub};
    stimuli = subj.stimuli;
    actions = subj.actions;
    outcomes = subj.outcomes;
    controllabilities = subj.controllability;

    schedule = 2 * (~subj.controllability(:, 1) & subj.isYoked(:, 1)) + subj.controllability(:, 1);

    for b = 1:nBlocks
        switch schedule(b)
            case 1
                hc = hc + 1;
                highControlStimuli(hc, :) = stimuli(b, :);
                highControlActions(hc, :) = actions(b, :);
                highControlOutcomes(hc, :) = outcomes(b, :);
            case 0
                lc = lc + 1;
                lowControlStimuli(lc, :) = stimuli(b, :);
                lowControlActions(lc, :) = actions(b, :);
                lowControlOutcomes(lc, :) = outcomes(b, :);
            case 2
                yc = yc + 1;
                yokedLowControlStimuli(yc, :) = stimuli(b, :);
                yokedLowControlActions(yc, :) = actions(b, :);
                yokedLowControlOutcomes(yc, :) = outcomes(b, :);
        end
    end
    
    % Initialize vectors
    maxWins = nTrials; % Maximum possible wins
    S = 4; % Number of stimuli
    
    % For high control
    shiftAfterLossCounts_HC = zeros(S, maxWins);
    totalConsecutiveWinsCounts_HC = zeros(S, maxWins);
    
    % For low control
    shiftAfterLossCounts_LC = zeros(S, maxWins);
    totalConsecutiveWinsCounts_LC = zeros(S, maxWins);
    
    % For yoked low control
    shiftAfterLossCounts_YC = zeros(S, maxWins);
    totalConsecutiveWinsCounts_YC = zeros(S, maxWins);
    
    % Process each type of control separately
    for controlType = {'HC', 'LC', 'YC'}
        ct = controlType{1};
        
        switch ct
            case 'HC'
                controlStimuli = highControlStimuli;
                controlActions = highControlActions;
                controlOutcomes = highControlOutcomes;
                shiftAfterLossCounts = shiftAfterLossCounts_HC;
                totalConsecutiveWinsCounts = totalConsecutiveWinsCounts_HC;
            case 'LC'
                controlStimuli = lowControlStimuli;
                controlActions = lowControlActions;
                controlOutcomes = lowControlOutcomes;
                shiftAfterLossCounts = shiftAfterLossCounts_LC;
                totalConsecutiveWinsCounts = totalConsecutiveWinsCounts_LC;
            case 'YC'
                controlStimuli = yokedLowControlStimuli;
                controlActions = yokedLowControlActions;
                controlOutcomes = yokedLowControlOutcomes;
                shiftAfterLossCounts = shiftAfterLossCounts_YC;
                totalConsecutiveWinsCounts = totalConsecutiveWinsCounts_YC;
        end
        
        for iS = 1:S
            for b = 1:size(controlStimuli, 1)
                consecutiveWins = 0;
                for t = 1:nTrials - 1
                    if controlStimuli(b, t) ~= iS
                        continue
                    end
                    s = controlStimuli(b, t);
                    o = controlOutcomes(b, t);
                    isWinState = mod(s, 2);
                    if (isWinState && o == 1) || (~isWinState && o == 0)
                        consecutiveWins = consecutiveWins + 1;
                    else
                        if consecutiveWins > 0
                            % Count the total number of consecutive wins
                            totalConsecutiveWinsCounts(iS, consecutiveWins) = totalConsecutiveWinsCounts(iS, consecutiveWins) + 1;                        
                            nextStateTrial = find(controlStimuli(b, t+1:end) == iS, 1, 'first') + t;
                            if ~isempty(nextStateTrial) && controlActions(b, t) ~= controlActions(b, nextStateTrial)
                                shiftAfterLossCounts(iS, consecutiveWins) = shiftAfterLossCounts(iS, consecutiveWins) + 1;
                            end
                        end
                        consecutiveWins = 0;
                    end
                end
            end
        end
        
        % Calculate probabilities
        probShiftAfterLoss = shiftAfterLossCounts ./ totalConsecutiveWinsCounts;
        
        % Handle division by zero (NaN values)
        probShiftAfterLoss(isnan(probShiftAfterLoss)) = 0;
        
        totalConsecutiveWinsCountsSum = sum(totalConsecutiveWinsCounts, 2);
        totalConsecutiveWinsCountsSum(totalConsecutiveWinsCountsSum == 0) = eps; % Avoid division by zero
        weightedProbShiftAfterLoss = probShiftAfterLoss .* totalConsecutiveWinsCounts ./ totalConsecutiveWinsCountsSum;
        
        % Store the results in the corresponding arrays
        switch ct
            case 'HC'
                weightedProbShiftAfterLoss_HC_perSubj = vertcat(weightedProbShiftAfterLoss_HC_perSubj, weightedProbShiftAfterLoss);
            case 'LC'
                weightedProbShiftAfterLoss_LC_perSubj = vertcat(weightedProbShiftAfterLoss_LC_perSubj, weightedProbShiftAfterLoss);
            case 'YC'
                weightedProbShiftAfterLoss_YC_perSubj = vertcat(weightedProbShiftAfterLoss_YC_perSubj, weightedProbShiftAfterLoss);
        end

        
    end
   
end

weightedProbShiftAfterLoss_HC = mean(vertcat(weightedProbShiftAfterLoss_HC_perSubj{:}), 1);
weightedProbShiftAfterLoss_LC = mean(vertcat(weightedProbShiftAfterLoss_LC_perSubj{:}), 1);
weightedProbShiftAfterLoss_YC = mean(vertcat(weightedProbShiftAfterLoss_YC_perSubj{:}), 1);

% Create a matrix with the three arrays for easier plotting
dataMatrix = [weightedProbShiftAfterLoss_HC; 
              weightedProbShiftAfterLoss_LC; 
              weightedProbShiftAfterLoss_YC]';

% Define the labels
labels = {'High Control', 'Low Control', 'Yoked Control'};

% Create a bar plot
figure;
bar(dataMatrix);
xlabel('Trial');
ylabel('Weighted Probability of Shift After Loss');
legend(labels, 'Location', 'Best');
title('Weighted Probability of Shift After Loss by Control Condition');
xlim([0, 10]);
grid on;




%% Accuracy by state and control type

controlTypes = {'High Control', 'Low Control', 'Yoked Control', 'All Data'};
controlData = {hc_d, lc_d, yc_d, data}; % Add data for all control types

meanGoBias = zeros(1, numel(controlTypes));
semGoBias = zeros(1, numel(controlTypes));

for iType = 1:numel(controlTypes)
    currentData = controlData{iType};
    nSub = numel(currentData);
    
    goBiasBySub = zeros(1, nSub);

    for iSub = 1:nSub
        subj = currentData{iSub};

        actions = subj.actions;
        states = subj.stimuli;
        goActions = actions == 1;
        noGoActions = actions == 2;
        accuracyGoToWin = mean(goActions(states == 1), 'all');
        accuracyNoGoToWin = mean(noGoActions(states == 3), 'all');
        goBiasBySub(iSub) = accuracyGoToWin - accuracyNoGoToWin;
    end

    meanGoBias(iType) = mean(goBiasBySub, 'omitnan');
    semGoBias(iType) = nanstd(goBiasBySub, [], 2) / sqrt(nSub);
end

% Plotting
figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
bar(meanGoBias, 'FaceColor', 'flat');
hold on;

% Add error bars
errorbar(1:numel(controlTypes), meanGoBias, semGoBias, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add labels and title
set(gca, 'XTickLabel', controlTypes);
xlabel('Control Type');
ylabel('Go Bias (P(Go|Win) - P(NoGo|Win))');
title('Go Bias by Control Type (Acc(GoToWin) - Acc(NoGoToWin))');
ylim([-1, 1]);
grid on;
hold off;


%% Best fitting for which model?
% This section should plot the average learning curve for hc, lc and yc for
% all subjects which are most fitting to a certain model.

% First, retrieve a list of subjects for each model.
% Loop over each sublist of subjects for each model, and filter the data.
% Plot the learning curves. 
% Plot other patterns if necessary.

% Initialize a variable to store the most recently modified file
mostRecentFile = '';
mostRecentDate = 0;
files = dir(dirs.hbi);
% Loop through each file to find the most recently modified one
for k = 1:length(files)
    % Check if the filename matches the pattern for model comparison files
    if contains(files(k).name, 'hbi_mod_') && contains(files(k).name, '_allData.mat')
        % Get the last modified date of the file
        fileDate = files(k).datenum;
        % Update the most recent file if this one is newer
        if fileDate > mostRecentDate
            mostRecentDate = fileDate;
            mostRecentFile = files(k).name;
        end
    end
end

% Display the most recently modified file
if ~isempty(mostRecentFile)
    fprintf('The most recently modified file is: \n%s\n', mostRecentFile);
else
    fprintf('No matching files found.\n');
end

f_hbi = load(mostRecentFile);
cbm = f_hbi.cbm;
nSub = size(cbm.output.parameters{1}, 1);
nMod = size(cbm.output.parameters, 1);
responsibility = cbm.output.responsibility;

subResp = nan(nSub, 1);
for iSub = 1:nSub
    subMax = max(responsibility(iSub, :));
    subResp(iSub) = find(responsibility(iSub, :) == subMax);
end
modCounter = 0;
modResp = cell(nMod, 1);

for iMod = 1:nMod
    modResp{iMod} = data(subResp == iMod);
end

% Groups of mods to be combined
group1Mods = [14];
group2Mods = [12];

% Combine data for groups
combinedGroup1 = concatenateGroupData(group1Mods, modResp);
combinedGroup2 = concatenateGroupData(group2Mods, modResp);

% Plot data for both groups
plotData(combinedGroup1, num2str(group1Mods, "M%02d_"));
plotData(combinedGroup2, num2str(group2Mods, "M%02d_"));

% How likely am I going to stay if I performed action Go/NoGo and got
% rewarded or punished? How does that differ between high and low control?
% Initialize P_stay to hold counts of stay actions and total trials for each condition
% Define data types

plotStayAnalysisForGroups(combinedGroup1, combinedGroup2, 'Group 1 (M14)', 'Group 2 (M12)');







%% ----------------------------------------------------------------------- %%
% Function to plot data using plotParticipantCurves
function plotData(combinedData, groupName)
    % Extract control type data
    [hc_d, lc_d, yc_d] = extractControlTypeData(combinedData);

    % Create a new figure and title it
    figure;
    sgtitle(['Learning curves for ' groupName]);

    % Subplot for hc_d
    subplot(1, 3, 1); % 1 row, 3 columns, 1st subplot
    plotParticipantCurves(hc_d, gcf);
    title('HC_D');

    % Subplot for lc_d
    subplot(1, 3, 2); % 1 row, 3 columns, 2nd subplot
    plotParticipantCurves(lc_d, gcf);
    title('LC_D');

    % Subplot for yc_d
    subplot(1, 3, 3); % 1 row, 3 columns, 3rd subplot
    plotParticipantCurves(yc_d, gcf);
    title('YC_D');
end

% Function to plot data for stay analysis for both groups
function plotStayAnalysisForGroups(combinedGroup1, combinedGroup2, groupName1, groupName2)
    % Create a new figure for Stay Analysis
    figure;
    sgtitle('Stay Analysis for Both Groups');

    % Subplot for stay analysis of the first group
    subplot(1, 2, 1); % 1 row, 2 columns, 1st subplot
    plotStayAnalysis(combinedGroup1, gcf);
    title(['Mean Probability of Staying - ' groupName1]);

    % Subplot for stay analysis of the second group
    subplot(1, 2, 2); % 1 row, 2 columns, 2nd subplot
    plotStayAnalysis(combinedGroup2, gcf);
    title(['Mean Probability of Staying - ' groupName2]);
end

% Function to concatenate data from specified groups
function combinedData = concatenateGroupData(groupMods, modResp)
    combinedData = horzcat(modResp{groupMods}); % Concatenate the cell arrays for specified mods
end

% Function to calculate the mean response time for a given cell array
function meanResponseTime = calculateMeanResponseTime(responseTimes)
    meanSize = floor(mean(cellfun(@length, responseTimes)));
    responseTimeMatrix = nan(numel(responseTimes), meanSize);

    for i = 1:numel(responseTimes)
        currentArray = responseTimes{i};
        if length(currentArray) >= meanSize
            responseTimeMatrix(i, :) = currentArray(1:meanSize);
        end
    end

    % Calculate the mean response time, ignoring NaNs
    meanResponseTime = nanmean(responseTimeMatrix, 1);
end