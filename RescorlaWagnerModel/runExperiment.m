function [blockInfo, HCprobGoMatrix, LCprobGoMatrix] = runExperiment(epsilon, beta, rho, numTrialsInBlock, numBlocks, rewardProb, controllProb)

    % Matrix to store info about blocks for later analysis
    blockInfo = zeros(numBlocks, 3);

    % Controllability array
    numHCBlocks = numBlocks / 2;
    controllabilityArray = zeros(numBlocks, 1);
    controllabilityArray(1:numHCBlocks) = 1;
    controllabilityArray = controllabilityArray(randperm(length(controllabilityArray)));
    
    conditions = [1, 2, 3, 4];
    repetitions = numTrialsInBlock / 4;

    % Model initialization
    m = Model(epsilon, rho, beta);

    currentTrial = 0;
    HCprobGoMatrix = cell(numHCBlocks, 4);
    LCprobGoMatrix = cell(numBlocks - numHCBlocks, 4); % Should be numLCBlocks, if these differ
    highControlCount = 0;
    lowControlCount = 0;

    for block = 1:numBlocks
        isHighControl = controllabilityArray(block);
        blockInfo(block, :) = [currentTrial + 1, currentTrial + numTrialsInBlock, isHighControl];
        m = m.resetP();
        m = m.resetQ();

        stateArray = repelem(conditions, repetitions);
        stateArray = stateArray(randperm(length(stateArray)));
        env = TrialEnvironment(rewardProb, controllProb, stateArray);

        if isHighControl
            highControlCount = highControlCount + 1;
        else
            lowControlCount = lowControlCount + 1;
        end

        for trial = 1:numTrialsInBlock
            [state, ~, env] = env.presentTrial();

            m = m.calcProbs(state);
            if isHighControl
                HCprobGoMatrix = updateProbGoMatrix(HCprobGoMatrix, highControlCount, state, m);
            else
                LCprobGoMatrix = updateProbGoMatrix(LCprobGoMatrix, lowControlCount, state, m);
            end

            action = m.returnAction(state);
            reward = env.getReward(state, action, isHighControl);
            m = m.updateModel(reward, state, action);
        end
    end
end

function matrix = updateProbGoMatrix(matrix, blockCount, state, model)
    if isempty(matrix{blockCount, state})
        matrix{blockCount, state} = model.returnGoProb(state);
    else
        matrix{blockCount, state}(end+1) = model.returnGoProb(state);
    end
end
