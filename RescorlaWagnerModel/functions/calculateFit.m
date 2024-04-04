function negLogLikelihood = calculateFit(model, simulatedActions, stimulusSequence, outcomeVec, successVec, correctActions)
    numBlocks = size(simulatedActions, 1);
    numTrials = size(simulatedActions, 2);
    negLogLikelihood = 0;

    for block = 1:numBlocks
        m = model;
        for trial = 1:numTrials
            state = stimulusSequence{block}(trial);
            simulatedAction = simulatedActions(block, trial);

            m = m.calcProbs(state);
            actionProb = m.P(state, simulatedAction);
            negLogLikelihood = negLogLikelihood - log(actionProb);

            action = m.returnAction(state); % Sample an action based on the model's predictions
            correctAction = correctActions(block, trial) == action;
            outcome = outcomeVec{block}(trial);
            success = successVec{block}(trial);
            reward = determineReward(outcome, success, state, correctAction);
            m = m.updateModel(reward, state, action);
        end
    end
end