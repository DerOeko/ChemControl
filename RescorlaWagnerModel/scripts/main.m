%% Main script of Rescorla Wagner model implementation
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
% Number of Trials per Block	| T	| The number of trials 
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
cfg = Config();

B = cfg.B; % Number of blocks
T = cfg.T; % Number of trials in blocks
R = cfg.R; % Number of Runs
model_names = cfg.model_names;

fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]); % Figure for non-Omega models
fig2 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
fig3 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]); % Figure for Omega models
fig4 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

numOmega = cfg.numOmegas; % How many omega models are there in total?
meanOmegas = zeros(B * T, numOmega); % Array of mean small omegas of each model
omegaCount = 0; % Counter for Omega Models (for plotting
nonOmegaCount = 0;

for i = 1:length(cfg.models)
    model = cfg.models{i};
    model_name = model_names{i};

    % ToDo fix average Experiment function
    [HCmeans, LCmeans, averageOmegas, cProbs] = averageExperiment(model, cfg);

    % Store the averageOmegas in the appropriate column of averageOmegasList
    if contains(model_name, "Fixed Omega Model", "IgnoreCase", true)
        meanOmegas(:, 1) = averageOmegas;
    elseif contains(model_name, "Dynamic Omega Model with Qoppa", "IgnoreCase", true)
        meanOmegas(:, 2) = averageOmegas;
    elseif contains(model_name, "Dynamic Omega Model with Omega", "IgnoreCase", true)
        meanOmegas(:, 3) = averageOmegas;
    % elseif contains(model_name, "Dynamic Omega Model with Associability", "IgnoreCase", true)
    %     averageOmegasList(:, 4) = averageOmegas;
    end

    if ~contains(model_name, "Omega", "IgnoreCase", true)
        nonOmegaCount = nonOmegaCount + 1; % Increment counter for non-Omega models
        fig = figure(fig1); % Focus on figure for non-Omega models
        subplot(2, 4, nonOmegaCount); % Plot in the first row for high control
        plotLearningCurves(HCmeans, model_name, true, fig);

        subplot(2, 4, nonOmegaCount + 4); % Plot in the second row for low control
        plotLearningCurves(LCmeans, model_name, false, fig);

    else
        omegaCount = omegaCount + 1; % Increment counter for Omega models
        fig = figure(fig3); % Focus on figure for Omega models
        subplot(2, numOmega, omegaCount); % Dynamic subplot placement for Omega models
        plotLearningCurves(HCmeans, model_name, true, fig);

        subplot(2, numOmega, omegaCount + numOmega);% Continue in the next row for low control
        plotLearningCurves(LCmeans, model_name, false, fig);
    end

    figure(fig2);
    subplot(2, ceil(size(model_names, 2)/2), i);
    plotBoxplots(HCmeans, LCmeans, model_name, fig2)

end

figure(fig1);
sgtitle('Comparison of Different RM Models');
legend('G2W', 'G2A', 'NG2W', 'NG2A', 'Location', 'best');

figure(fig2);
ylabel('Proportion of Bias');
xlabel('Condition');

figure(fig3);
sgtitle('Comparison of Omega RM Models');
legend('G2W', 'G2A', 'NG2W', 'NG2A', 'Location', 'best');

fig = figure(fig4);
plotOmegas(meanOmegas, cProbs, cfg, fig)
