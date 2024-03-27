%% Implementation of a Rescorla Wagner model with feedback sensitivity for Robot GoNoGo task
% Author: Samuel Nellessen, MCC

% Clear previous sessions for a clean start
close all;
clear all;

% Initialize hyperparameters
epsilon = 0.2;
beta = 3;
rho = 0.3;
numTrialsInBlock = 40;
numBlocks = 28;
rewardProb = 0.85;
controllProb = 0.8;
numHCBlocks = numBlocks / 2;
numLCBlocks = numBlocks - numHCBlocks;


[blockInfo, HCprobGoMatrix, LCprobGoMatrix] = runExperiment(epsilon, beta, rho, numTrialsInBlock, numBlocks, rewardProb, controllProb);

% First, you get all state == 1 in a block. the row should be block number
% THe column should be the condition repetition, with each cell being the probability of choosing "Go" in that condition.
% Then, you average over all the columns, to get the consolidated graph. 
% Then you plot this.
% Repeat for all conditions
% block, state, say block 1, and state 1. Then block 1, state 1 is the
% probabilities of choosing Go for each trial where it was state 1 in block
% 1. 
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


lineStyles = {'-', '--', ':', '-.'}; % Example line styles
markers = {'o', '+', '*', 'x'}; % Example markers
colors = {'r', 'g', 'b', 'c'}; % Example colors (red, green, blue, purple)
figure;
subplot(2,1,1);
hold on; % Allows multiple plots on the same figure

for state = 1:4
    plotStyle = [colors{state} lineStyles{state} markers{state}]; % Combine color, line style, and marker
    plot(1:10, HCoccurrenceMeans(:, state)', plotStyle, 'LineWidth', 2); % Plotting mean probabilities for each state
end

xlabel('State Repetitions');
ylabel('P(Go response | state)');
ylim([0.0, 1.0])
yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
legend('GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss', 'Location', 'best');
title('Mean Probability of Choosing "Go" Across State Repetitions In High Control Trials');
hold off;

subplot(2,1,2);
hold on;

for state = 1:4
    plotStyle = [colors{state} lineStyles{state} markers{state}]; % Combine color, line style, and marker
    plot(1:10, LCoccurrenceMeans(:, state)', plotStyle, 'LineWidth', 2); % Plotting mean probabilities for each state
end

xlabel('State Repetitions');
ylabel('P(Go response | state)');
ylim([0.0, 1.0])
yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
legend('GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss', 'Location', 'best');
title('Mean Probability of Choosing "Go" Across State Repetitions in Low Control Trials');
hold off;


% Calculate Pavlovian bias
GoToWin_GoToAvoid_HC = HCoccurrenceMeans(:, 1) - HCoccurrenceMeans(:, 2);
% Reasoning behind the 1-mean here:
% Pavlovian bias (also called motivational bias) denotes that phenomenon that 
% - reward-related cues (eliciting reward anticipation) invigorate action 
% (lead to more active "Go" responses and speed up these Go responses) 
% - punishment-related cues (eliciting punishment anticipation) suppress action 
% (lead to less "Go" / more "NoGo" responses and slow down Go responses).
% Taken from cognitive atlas.
% With this idea, I thought that we want to show the bias in seeking
% rewards generally, as well as the bias in avoiding losses generally.
% The former is the mean probability across states and acrsss blocks to
% choose "Go" in a win condition, minus the probability of choosing "Go" in
% a loss condition. Ideally, they should be the same (i.e., an ideal learner
% should show a Go response in either case, because it is the correct
% action). Intuitively, this shows how much more likely we are seeking
% rewards (we react stronger to rewards). 
% Conversely, we want to know how much more likely we are avoiding losses.
% So, when "NoGo" is the correct action, how much higher is the probability
% to "NoGo" when the condition is loss than when the condition is win
% (again, for an ideal learner, this would be 0). 
% To compute this, we need the probability of choosing "NoGo" when the
% setting is a loss setting and the correct action is "NoGo", and when the
% setting is win and the correct action is "NoGo". 
NoGoToAvoid_NoGoToWin_HC = (1-HCoccurrenceMeans(:, 4)) - (1-HCoccurrenceMeans(:, 3));

% Calculate Pavlovian bias
GoToWin_GoToAvoid_LC = LCoccurrenceMeans(:, 1) - LCoccurrenceMeans(:, 2);
NoGoToAvoid_NoGoToWin_LC = (1-LCoccurrenceMeans(:, 4)) - (1-LCoccurrenceMeans(:, 3));

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

% 
% % Plotting
% figure;
% plot(1:length(meanProbabilities), meanProbabilities, '-o', 'LineWidth', 2);
% title(sprintf('Mean Probability of Choosing "Go" for State %d Across Blocks', chosenState));
% xlabel('Occurrence Number');
% ylabel('Mean Probability');
% grid on;
% proportionCorrectHigh = zeros(numBlocks, 1);
% proportionCorrectLow = zeros(numBlocks, 1);
% rewardInHighBlocks = zeros(numBlocks,1);
% rewardInLowBlocks = zeros(numBlocks, 1);
% 
% countHigh = 0;
% countLow = 0;
% 
% for block = 1:numBlocks
%     startIndex = (block-1) * numTrialsInBlock +1;
%     endIndex = (block)*numTrialsInBlock;
%     isHighControl = blockInfo(block,3);
% 
%     if isHighControl
%         countHigh = countHigh + 1;
%         proportionCorrectHigh(countHigh) = sum(correctAnswers(startIndex:endIndex) == 1) / numTrialsInBlock;
%         rewardInHighBlocks(countHigh) = rewardInBlocks(block);
%     else
%         countLow = countLow + 1;
%         proportionCorrectLow(countLow) = sum(correctAnswers(startIndex:endIndex) == 1) / numTrialsInBlock;
%         rewardInLowBlocks(countLow) = rewardInBlocks(block);
%     end
% end
% 
% % Trim the arrays to remove unused elements (if any)
% proportionCorrectHigh = proportionCorrectHigh(1:countHigh);
% proportionCorrectLow = proportionCorrectLow(1:countLow);
% rewardInHighBlocks = rewardInHighBlocks(1:countHigh);
% rewardInLowBlocks = rewardInLowBlocks(1:countLow);
% 
% % Plot for High Control Blocks
% figure;
% bar(proportionCorrectHigh, 'FaceColor', 'blue');
% title('Proportion of Correct Answers in High Control Blocks');
% xlabel('Block Number');
% ylabel('Proportion of Correct Answers');
% ylim([0, 1]);
% 
% % Plot for Low Control Blocks
% figure;
% bar(proportionCorrectLow, 'FaceColor', 'red');
% title('Proportion of Correct Answers in Low Control Blocks');
% xlabel('Block Number');
% ylabel('Proportion of Correct Answers');
% ylim([0, 1]);
% 
% % Plot for Rewards in High Control Blocks
% figure;
% plot(1:length(rewardInHighBlocks), rewardInHighBlocks, 'LineWidth', 3, 'Color', 'blue');
% title('Reward in High Control Blocks');
% xlabel('Block Number');
% ylabel('Total Reward');
% 
% % Plot for Rewards in Low Control Blocks
% figure;
% plot(1:length(rewardInLowBlocks), rewardInLowBlocks, 'LineWidth', 3, 'Color', 'red');
% title('Reward in Low Control Blocks');
% xlabel('Block Number');
% ylabel('Total Reward');
% 
% 
% % Smooth the correct answers data using a moving average and plot it
% smoothedCorrectAnswers = movmean(correctAnswers, 10);
% trialNumbers = 1:length(correctAnswers);
% figure;
% plot(trialNumbers, smoothedCorrectAnswers, 'LineWidth', 2);
% title('Running Average of Correct Answers');
% xlabel('Trial');
% ylabel('Proportion of Correct Answers');
% ylim([-0.01, 1.01]);
% 
% % Plot markers for correct and incorrect answers
% hold on;
% plot(trialNumbers(correctAnswers == 1), 0.99*ones(size(trialNumbers(correctAnswers == 1))), 'p', 'MarkerSize', 5, 'MarkerEdgeColor', 'black', 'LineStyle', 'none');
% plot(trialNumbers(correctAnswers == 0), 0.01*ones(size(trialNumbers(correctAnswers == 0))), 'p', 'MarkerSize', 5, 'MarkerEdgeColor', 'black', 'LineStyle', 'none');
% 
% % Add shaded areas to indicate high and low controllability blocks
% highControlColor = [0.8, 0.8, 1]; % Light blue for high control
% lowControlColor = [1, 0.8, 0.8]; % Light red for low control
% 
% for block = 1:numBlocks
%     startTrial = blockInfo(block, 1);
%     endTrial = blockInfo(block, 2);
%     isHighControl = blockInfo(block, 3);
% 
%     x = [startTrial, endTrial, endTrial, startTrial];
%     y2 = [min(cumulativeRewards), min(cumulativeRewards), max(cumulativeRewards), max(cumulativeRewards)];
% 
%     if isHighControl
%         color = highControlColor;
%     else
%         color = lowControlColor;
%     end    
%     fill(x, y2, color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% end
% 
% hold off; % End of the plotting section
% 
% % Calculate average accuracy for high control blocks
% averageAccuracyHigh = mean(proportionCorrectHigh);
% 
% % Calculate average accuracy for low control blocks
% averageAccuracyLow = mean(proportionCorrectLow);
% 
% % Print out the average accuracies
% fprintf('Average Accuracy for High Control Blocks: %.2f%%\n', averageAccuracyHigh * 100);
% fprintf('Average Accuracy for Low Control Blocks: %.2f%%\n', averageAccuracyLow * 100);

