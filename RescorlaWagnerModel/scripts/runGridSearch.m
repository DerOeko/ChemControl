function [accuracyMatrix, epsilonFlat, betaFlat, rhoFlat, accuracyFlat] = runGridSearch(epsilons, betas, rhos, numTrialsInBlock, numBlocks, rewardProb, controllProb)
    % Validate or set default parameter values for epsilons, betas, rhos, etc., if not passed

    arguments
        % Parameter ranges
        epsilons = 0.01:0.02:0.8;
        betas = 1:1:15;
        rhos = 0.01:1:15;

        % Set the number of trials per block and the number of blocks
        numTrialsInBlock = 40;
        numBlocks = 16;
    
        % Initialize the environment with success probabilities for the two states
        rewardProb = 0.85;
        controllProb = 0.8;
    end
    
    % Controllability array
    numHCBlocks = numBlocks/2;
    numLCBlocks = numBlocks/2;
    
    conditions = [1,2,3,4];
    repetitions = numTrialsInBlock/4;
    

    results = zeros(length(epsilons), length(betas), length(rhos));
    accuracyMatrix = zeros(length(epsilons), length(betas), length(rhos));
    currentTrial = 0;
    for i = 1:length(epsilons)
        for j = 1:length(betas)
            for k = 1:length(rhos)
                epsilon = epsilons(i);
                beta = betas(j);
                rho = rhos(k);
                
                controllabilityArray = zeros(numBlocks);
                controllabilityArray(1:numHCBlocks) = 1;
                controllabilityArray = controllabilityArray(randperm(length(controllabilityArray)));
                % Initialize indices for tracking trials and total reward accumulated
                trialIndex = 0;
                % Variables to store correct actions and cumulative rewards for analysis
                correctAnswers = zeros(numTrialsInBlock*numBlocks, 1);
    
                highControlAccuracies = []; % Collect accuracies for high control blocks
                m = Model(epsilon, rho, beta);
    
                % Begin experiment by iterating through each block
                for block = 1:numBlocks
                    isHighControl = controllabilityArray(block);
    
                    correctAnswers = zeros(numTrialsInBlock, 1); % Reset for each block
    
                    m = m.resetP();
                    m = m.resetQ();
                
                    stateArray = repelem(conditions, repetitions);
                    stateArray = stateArray(randperm(length(stateArray)));
                    env = TrialEnvironment(rewardProb, controllProb, stateArray);
    
    
                    % Iterate through trials within the current block
                    for trial = 1:numTrialsInBlock
                        trialIndex = trialIndex + 1;
                        [state, correctAction, env] = env.presentTrial();
                
                        m = m.calcProbs(state);
                        % Present the current trial and get the state and correct action
    
    
                        action = m.returnAction(state);        
                        reward = env.getReward(state, action, isHighControl);
                        m = m.updateModel(reward, state, action);
    
                        % Record whether the chosen action was correct and update the total reward
                        correctAnswers(trialIndex) = (action == correctAction);
                    end
                    if isHighControl
                        blockAccuracy = sum(correctAnswers) / numTrialsInBlock; % Calculate for current block only
                        highControlAccuracies = [highControlAccuracies, blockAccuracy]; % Append to array
                    end
                end
                accuracyMatrix(i, j, k) = mean(highControlAccuracies);
            end
        end
    end

    % Flatten the parameter grids and accuracies
    [epsilonGrid, betaGrid, rhoGrid] = ndgrid(epsilons, betas, rhos);
    epsilonFlat = epsilonGrid(:);
    betaFlat = betaGrid(:);
    rhoFlat = rhoGrid(:);
    accuracyFlat = accuracyMatrix(:);
end