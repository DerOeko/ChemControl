% Clear previous sessions for a clean start
close all;
clear all;

% Initialize learning rate (epsilon) for Q-value updates
epsilon = 0.1;

% Initialize softmax temperature (beta) for the softmax action selection
beta = 5;

% Feedback sensitivity (rho) scales the difference between received and estimated reward
rho = 10;

% Set the number of trials per block and the number of blocks
numTrialsInBlock = 40;
numBlocks = 16;

% Initialize action values (Q-values) for 4 state-action pairs and the action probability matrix
Q = zeros(4,2);
p = zeros(4,2);

% Initialize the environment with success probabilities for the two states
env = TrialEnvironment(0.85, 0.8);

% Variables to store correct actions and cumulative rewards for analysis
correctAnswers = zeros(numTrialsInBlock*numBlocks, 1);
cumulativeRewards = zeros(numTrialsInBlock*numBlocks, 1);

% Initialize indices for tracking trials and total reward accumulated
trialIndex = 0;
totalReward = 0;

% Probability of being in a high control block
highControlProb = 0.5;

% Default controllability setting for blocks
controlString = "high";

% Matrix to store info about blocks for later analysis
blockInfo = zeros(numBlocks, 3);

% Track the current trial number globally
currentTrial = 0;

% Begin experiment by iterating through each block
for block = 1:numBlocks
    % Determine the controllability condition for the current block randomly
    isHighControl = rand < highControlProb;

    % Set the control string based on the condition
    if isHighControl
        controlString = "high";
    else
        controlString = "low";
    end
    % Store the block information including start and end trial numbers, and controllability condition
    blockInfo(block, :) = [currentTrial + 1, currentTrial + numTrialsInBlock, isHighControl];

    % Print block information
    fprintf('\n================ Block %d: Controllability %s ==========\n', block, controlString);
    
    % Iterate through trials within the current block
    for trial = 1:numTrialsInBlock
        trialIndex = trialIndex + 1;
        
        % Present the current trial and get the state and correct action
        [state, correctAction] = env.presentTrial();
        
        % Display trial state and correct action
        fprintf('State: %d\n', state);
        fprintf('Correct Action: %d\n', correctAction);
        
        % Calculate action probabilities using the softmax function
        p(state, :) = beta_softmax(Q(state, :), beta);
        
        % Select an action based on calculated probabilities
        action = randsample([1,2], 1, true, p(state, :));
        
        % Display chosen action
        fprintf('Chosen Action: %d\n', action);
        
        % Get the reward for the chosen action in the current state and controllability condition
        reward = env.getReward(state, action, isHighControl);
        
        % Update the Q-value for the chosen action using the learning rate and feedback sensitivity
        Q(state, action) = Q(state, action) + epsilon * (rho * reward - Q(state, action)); 
        
        % Record whether the chosen action was correct and update the total reward
        correctAnswers(trialIndex) = (action == correctAction);
        totalReward = totalReward + reward;
        cumulativeRewards(trialIndex) = totalReward;
        currentTrial = currentTrial + 1;
    end
end

% Calculate the proportion of correct answers per block for analysis
proportionCorrect = zeros(numBlocks, 1);
for block = 1:numBlocks
    startIndex = (block - 1) * numTrialsInBlock + 1;
    endIndex = block * numTrialsInBlock;
    proportionCorrect(block) = sum(correctAnswers(startIndex:endIndex) == 1) / numTrialsInBlock;
end

% Plot the proportion of correct answers by block
figure;
bar(proportionCorrect, 'FaceColor', 'blue');
title('Proportion of Correct Answers by Block');
xlabel('Block Number');
ylabel('Proportion of Correct Answers');
ylim([0, 1]);

% Smooth the correct answers data using a moving average and plot it
smoothedCorrectAnswers = movmean(correctAnswers, 5);
trialNumbers = 1:length(correctAnswers);
figure;
plot(trialNumbers, smoothedCorrectAnswers, 'LineWidth', 2);
title('Running Average of Correct Answers');
xlabel('Trial');
ylabel('Proportion of Correct Answers');
ylim([-0.01, 1.01]);

% Plot markers for correct and incorrect answers
hold on;
plot(trialNumbers(correctAnswers == 1), 0.99*ones(size(trialNumbers(correctAnswers == 1))), 'p', 'MarkerSize', 5, 'MarkerEdgeColor', 'black', 'LineStyle', 'none');
plot(trialNumbers(correctAnswers == 0), 0.01*ones(size(trialNumbers(correctAnswers == 0))), 'p', 'MarkerSize', 5, 'MarkerEdgeColor', 'black', 'LineStyle', 'none');

% Add shaded areas to indicate high and low controllability blocks
highControlColor = [0.8, 0.8, 1]; % Light blue for high control
lowControlColor = [1, 0.8, 0.8]; % Light red for low control

for block = 1:numBlocks
    startTrial = blockInfo(block, 1);
    endTrial = blockInfo(block, 2);
    isHighControl = blockInfo(block, 3);

    x = [startTrial, endTrial, endTrial, startTrial];
    y2 = [min(cumulativeRewards), min(cumulativeRewards), max(cumulativeRewards), max(cumulativeRewards)];

    if isHighControl
        color = highControlColor;
    else
        color = lowControlColor;
    end    
    fill(x, y2, color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

hold off; % End of the plotting section
