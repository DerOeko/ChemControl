%% Implementation of a Rescorla Wagner model with feedback sensitivity for Robot GoNoGo task
% Author: Samuel Nellessen, MCC

% Clear previous sessions for a clean start
close all;
clear all;

lineStyles = {'-', '--', ':', '-.'}; % Example line styles
markers = {'o', '+', '*', 'x'}; % Example markers
colors = {'r', 'g', 'b', 'c'}; % Example colors (red, green, blue, purple)

% Initialize hyperparameters
epsilon = 0.2;
beta = 1;
rho = 1;
numTrialsInBlock = 40;
numBlocks = 28;
rewardProb = 0.85;
controllProb = 0.8;
numHCBlocks = numBlocks / 2;
numLCBlocks = numBlocks - numHCBlocks;
numIterations = 1;

% Initializing accumulators for High Control (HC) and Low Control (LC)
accumHCMeans = zeros(10, 4);  
accumLCMeans = zeros(10, 4);

for iter= 1:numIterations
    [blockInfo, HCprobGoMatrix, LCprobGoMatrix] = runExperiment(epsilon, beta, rho, numTrialsInBlock, numBlocks, rewardProb, controllProb);
    
    HCoccurrenceMeans = NaN(10, 4);
    for state = 1:4
        for occurrence = 1:10
            HCoccurrenceProbs = zeros(numHCBlocks, 1);
            for block = 1:numHCBlocks
                HCoccurrenceProbs(block) = HCprobGoMatrix{block, state}(occurrence);
            end
            HCoccurrenceMeans(occurrence, state) = mean(HCoccurrenceProbs);
        end
    end
    
    LCoccurrenceMeans = NaN(10, 4);
    for state = 1:4
        for occurrence = 1:10
            LCoccurrenceProbs = zeros(numLCBlocks, 1);
            for block = 1:numLCBlocks
                LCoccurrenceProbs(block) = LCprobGoMatrix{block, state}(occurrence);
            end
            LCoccurrenceMeans(occurrence, state) = mean(LCoccurrenceProbs);
        end
    end

    accumHCMeans = accumHCMeans + HCoccurrenceMeans;  
    accumLCMeans = accumLCMeans + LCoccurrenceMeans;
end

meanHCMeans = accumHCMeans / numIterations;
meanLCMeans = accumLCMeans / numIterations;

figure;
subplot(2,1,1);
hold on; % Allows multiple plots on the same figure

for state = 1:4
    plotStyle = [colors{state} lineStyles{state} markers{state}]; % Combine color, line style, and marker
    plot(1:10, meanHCMeans(:, state)', plotStyle, 'LineWidth', 2); % Plotting mean probabilities for each state
end

xlabel('State Repetitions');
ylabel('P(Go response | state)');
ylim([0.0, 1.0])
yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
legend('GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss', 'Location', 'best');
title('Mean Probability of Choosing "Go" Across State Repetitions In High Control Trials Across 50 Experiments');
hold off;

subplot(2,1,2);
hold on;

for state = 1:4
    plotStyle = [colors{state} lineStyles{state} markers{state}]; % Combine color, line style, and marker
    plot(1:10, meanLCMeans(:, state)', plotStyle, 'LineWidth', 2); % Plotting mean probabilities for each state
end

xlabel('State Repetitions');
ylabel('P(Go response | state)');
ylim([0.0, 1.0])
yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
legend('GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss', 'Location', 'best');
title('Mean Probability of Choosing "Go" Across State Repetitions in Low Control Trials Across 50 Experiments');
hold off;


% Calculate Pavlovian bias
GoToWin_GoToAvoid_HC = meanHCMeans(:, 1) - meanHCMeans(:, 2);
NoGoToAvoid_NoGoToWin_HC = (1-meanHCMeans(:, 4)) - (1-meanHCMeans(:, 3));

% Calculate Pavlovian bias
GoToWin_GoToAvoid_LC = meanLCMeans(:, 1) - meanLCMeans(:, 2);
NoGoToAvoid_NoGoToWin_LC = (1-meanLCMeans(:, 4)) - (1-meanLCMeans(:, 3));

% Single matrix for boxplot
biases = [GoToWin_GoToAvoid_HC, GoToWin_GoToAvoid_LC, NoGoToAvoid_NoGoToWin_HC, NoGoToAvoid_NoGoToWin_LC];
figure;
boxplot(biases, 'Labels', {'GoToWin-GotoAvoid in HC', 'GoToWin-GotoAvoid in LC' 'NoGoToAvoid-NoGoToWin in HC', 'NoGoToAvoid-NoGoToWin in LC'});
hold on;

numDataPoints = size(biases, 1);
% Create scatter plot for each group
for i = 1:4
    scatter(repelem(i, numDataPoints), biases(:, i), 'd', 'filled');
end

hold off;
title('Pavlovian Bias in High and Low Control Trials');
ylabel('Proportion of Bias');
xlabel('Condition');


