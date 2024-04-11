%% Implementation of a Rescorla Wagner model with feedback sensitivity for Robot GoNoGo task
% Author: Samuel Nellessen, MCC

%% Parameter dictionary
% Only for explanation of the used parameters, mostly for myself.

% Parameter name | In code | Conceptually
% Learning Rate	| epsilon	 | The rate at which the model updates its 
% Q-values based on the difference between the expected and actual reward. 
% A higher epsilon leads to faster learning.
% Feedback Sensitivity	| rho	| The scaling factor that determines how 
% much the experienced reward influences the Q-value update. A higher rho 
% means the model is more sensitive to the feedback.
% Softmax Temperature	| beta	|The temperature parameter that controls 
% the exploration-exploitation trade-off in the softmax action selection. 
% A higher beta leads to more deterministic, exploitative behavior.
% Number of Trials per Block	| numTrialsInBlock	| The number of trials 
% in each experimental block.
% Number of Blocks	| numBlocks	| The total number of blocks in the 
% experiment.
% Reward Probability	| rewardProb	| The probability that the correct 
% action leads to a good outcome (reward) in the task.
% Controllability Probability	| controllProb	| The probability that 
% the outcome of the trial depends on the agent's action (i.e., the action 
% matters).
% Number of High Control Blocks |	numHCBlocks	 |The number of blocks 
% where the outcome depends on the agent's action (high controllability).
% Number of Low Control Blocks	| numLCBlocks	| The number of blocks 
% where the outcome does not depend on the agent's action (low controllability).
% Go Bias	| goBias | 	A bias added to the "Go" action value, which makes 
% the model more likely to choose the "Go" action.
% Motivational Bias Vector	| V	| A vector representing the motivational 
% bias for each state (e.g., reward-related states vs. punishment-related states).
% Motivational Bias Scaling Factor	| pi | A scaling factor that determines 
% the strength of the motivational bias.
% Action values | Q | Found in model definition. 4,2 matrix that updates
% the value of an action in a given state.
% Action weights | W | Found in model definition. 4,2 matrix, whose first
% column incorporates the go bias and motivational bias.
% States or Conditions | conditions/state | Each state the agent is going
% through. Can be Go2Win, Go2Avoid, NoGo2Win, NoGo2Avoid (in that order as
% well).
% Action propensities/probabilities | P | Probability of choosing each
% action in a given state. Each row sums to 1.
% Softmax function | betaSoftmax | Basically, looks at the ratio between Go
% and NoGo for a given state, and transforms the ratio into a proper
% probability distribution, i.e. strictly continuous between 0 and 1.


% Clear previous sessions for a clean start
close all;
clear all;
clc;

lineStyles = {'-', '--', ':', '-.'}; % Plotting line styles

% Initialize hyperparameters
epsilon = 0.31;
beta = 1;
rho = 0.5;
numTrialsInBlock = 80;
numBlocks = 8;
rewardProb = 0.85;
controllProb = 0.8;
numHCBlocks = numBlocks / 2;
numLCBlocks = numBlocks - numHCBlocks;
goBias = 0.3;
Qinit = [0.5*rho -0.5*rho 0.5*rho -0.5*rho; 0.5*rho -0.5*rho 0.5*rho -0.5*rho]';
V = [0.3 -0.3 0.3 -0.3];
Vinit = [0.3 -0.3 0.3 -0.3];
pi = 0.3;
numRuns = 100;

models = {
    Model(epsilon, rho, beta, Qinit),
    GoBiasModel(epsilon, rho, beta, Qinit, goBias),
    FixedMotivationalBiasModel(epsilon, rho, beta, Qinit, goBias, V, pi),
    DynamicMotivationalBiasModel(epsilon, rho, beta, Qinit, goBias, Vinit, pi)
};

model_names = {'Generic Model', 'Go Bias Model', 'Fixed Motivational Bias Model', 'Dynamic Motivational Bias Model'};


%% Model definition
% choose either
% Model 1: Basic Model:
%model = Model2(epsilon, rho, beta);

% Model 2: Basic Model + Go Bias:
%model = GoBiasModel(epsilon, rho, beta, goBias);

% Model 3: Basic Model + Go Bias + Motivational Bias:
%model = FixedMotivationalBiasModel(epsilon, rho, beta, goBias, V, pi);

% Model 4: Basic Model + Go Bias + Dynamic Motivational Bias:
% model = DynamicMotivationalBiasModel(epsilon, rho, beta, goBias, Vinit, pi)
figure;
betas = {1 2 3 4 5};
for i = 1:5
    beta = betas{i};
    model = Model2(epsilon, rho, beta, Qinit);
    
    model_name = "Basic";
    [HCoccurrenceMeans, LCoccurrenceMeans] = averageExperiment(numRuns, model, numTrialsInBlock, numBlocks, rewardProb, controllProb);
    
    subplot(2,5,i)
    hold on; % Allows multiple plots on the same figure
    
    for state = 1:4
        plotStyle = [lineStyles{state} ]; % Combine color, line style, and marker
        plot(1:numTrialsInBlock/4, HCoccurrenceMeans(:, state)', 'LineWidth', 2); % Plotting mean probabilities for each state
    end
    
    xlabel('State Repetitions');
    xlim([1.0 numTrialsInBlock/4])
    ylabel('P(Go response | state)');
    ylim([0.0, 1.0])
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
    legend('G2W', 'G2A', 'NG2W', 'NG2A', 'Location', 'best');
    title(sprintf('%s \n with beta %i: \nMean P(Go|State) \n Across State Repetitions \nin High Control Trials', model_name, beta));
    grid on
    hold off;
    
    subplot(2, 5, i+5);
    hold on;
    
    for state = 1:4
        plotStyle = [ lineStyles{state}]; % Combine color, line style, and marker
        plot(1:numTrialsInBlock/4, LCoccurrenceMeans(:, state)', 'LineWidth', 2); % Plotting mean probabilities for each state
    end
    
    xlabel('State Repetitions');
    ylabel('P(Go response | state)');
    xlim([1.0 numTrialsInBlock/4])
    ylim([0.0, 1.0])
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
    title(sprintf('%s\n with beta %i: \nMean P(Go|State) \n Across State Repetitions \nin Low Control Trials', model_name, beta));
    grid on
    hold off;
end

%{
figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for i = 1:length(models)
    model = models{i};
    model_name = model_names{i};

    [HCoccurrenceMeans, LCoccurrenceMeans] = averageExperiment(numRuns, model, numTrialsInBlock, numBlocks, rewardProb, controllProb);
    disp("Done!")

    % Create the subplot for the current model
    subplot(2, 4, i);
    hold on; % Allows multiple plots on the same figure

    for state = 1:4
        plotStyle = [lineStyles{state} ]; % Combine color, line style, and marker
        plot(1:numTrialsInBlock/4, HCoccurrenceMeans(:, state)', 'LineWidth', 2); % Plotting mean probabilities for each state
    end

    xlabel('State Repetitions');
    xlim([1.0 numTrialsInBlock/4])
    ylabel('P(Go response | state)');
    ylim([0.0, 1.0])
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
    legend('G2W', 'G2A', 'NG2W', 'NG2A', 'Location', 'best');
    title(sprintf('%s: \nMean P(Go|State) \n Across State Repetitions \nin High Control Trials', model_name));
    grid on
    hold off;

    subplot(2, 4, i + 4);
    hold on;

    for state = 1:4
        plotStyle = [ lineStyles{state}]; % Combine color, line style, and marker
        plot(1:numTrialsInBlock/4, LCoccurrenceMeans(:, state)', 'LineWidth', 2); % Plotting mean probabilities for each state
    end

    xlabel('State Repetitions');
    ylabel('P(Go response | state)');
    xlim([1.0 numTrialsInBlock/4])
    ylim([0.0, 1.0])
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
    title(sprintf('%s: \nMean P(Go|State) \n Across State Repetitions \nin Low Control Trials', model_name));
    grid on
    hold off;
end


sgtitle('Comparison of Different RM Models');
legend('G2W', 'G2A', 'NG2W', 'NG2A', 'Location', 'best');

%}
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
NoGoToAvoid_NoGoToWin_HC = (HCoccurrenceMeans(:, 4)) - (HCoccurrenceMeans(:, 3));

% Calculate Pavlovian bias
GoToWin_GoToAvoid_LC = LCoccurrenceMeans(:, 1) - LCoccurrenceMeans(:, 2);
NoGoToAvoid_NoGoToWin_LC = (LCoccurrenceMeans(:, 4)) - (LCoccurrenceMeans(:, 3));

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

