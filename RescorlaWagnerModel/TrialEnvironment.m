classdef TrialEnvironment
    properties
        successProb = 0.85;
        outcomeProb = 0.8;
        % Assuming env.conditions = {'GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss'}
        % and actions are represented as [1, 2] for ['Go', 'NoGo']
        conditions = {1,2,3,4};
        stateArray;
        numTrialsInBlock;
        outcomeVec;
        successVec;
    end

    methods
        function obj = TrialEnvironment(successProb, outcomeProb, stateArray, numTrialsInBlock)
            % How likely is it that the correct action leads to a good
            % outcome?
            obj.successProb = successProb;
            % How likely is it that the action matters?
            obj.outcomeProb = outcomeProb;

            obj.stateArray = stateArray;
            nOutcomes = round(outcomeProb * numTrialsInBlock);
            obj.outcomeVec = [ones(1, nOutcomes) zeros(1, numTrialsInBlock-nOutcomes, 1)];
            obj.outcomeVec = obj.outcomeVec(randperm(numTrialsInBlock));
            nSuccesses = round(successProb * numTrialsInBlock);
            obj.successVec = [ones(1, nSuccesses) zeros(1, numTrialsInBlock-nSuccesses, 1)];
            obj.successVec = obj.successVec(randperm(numTrialsInBlock));
        end
        
        function [state, correctAction, obj] = presentTrial(obj)
            state = obj.stateArray(1);
            obj.stateArray(1) = [];
            correctAction = obj.getCorrectAction(state);
        end

        function correctAction = getCorrectAction(obj, state)
            if state <= 2
                correctAction = 1;
            else
                correctAction = 2;
            end
        end

        function reward = getReward(obj, trialIndex, state, action, isHighControl)

            outcomeMatters = obj.outcomeVec(trialIndex);
            successOutcome = obj.successVec(trialIndex);
            isLossState = ~mod(state, 2);
            correctAction = action == obj.getCorrectAction(state);
            fi = @(varargin) varargin{length(varargin) - varargin{1}};
        
            if isHighControl
                if outcomeMatters
                    if correctAction
                        if isLossState
                            reward = fi(successOutcome, 0, -10);
                        else
                            reward = fi(successOutcome, 10, 0);
                        end
                    else
                        if isLossState
                            reward = fi(~successOutcome, 0, -10);
                        else
                            reward = fi(~successOutcome, 10, 0);
                        end
                    end
                else
                    if isLossState
                        reward = fi(successOutcome, 0, -10);
                    else
                        reward = fi(successOutcome, 10, 0);
                    end
                end
            else
                if ~outcomeMatters
                    if correctAction
                        if isLossState
                            reward = fi(successOutcome, 0, -10);
                        else
                            reward = fi(successOutcome, 10, 0);
                        end
                    else
                        if isLossState
                            reward = fi(~successOutcome, 0, -10);
                        else
                            reward = fi(~successOutcome, 10, 0);
                        end
                    end
                else
                    if isLossState
                        reward = fi(successOutcome, 0, -10);
                    else
                        reward = fi(successOutcome, 10, 0);
                    end
                end
            end
            

                    reward = fi(correctAction && isLossState)
                    reward = fi(correctAction, (isLossState && successOutcome) * 0, ...
                                successOutcome * 10, ...
                                (isLossState && ~successOutcome) * 0, ...
                                ~successOutcome * 10);
                else
                    reward = fi(isLossState, (successOutcome * 0), (successOutcome * -10), ...
                                (~successOutcome * 0), (~successOutcome * -10));
                end
            else
                if ~outcomeMatters
                    reward = fi(correctAction, (isLossState && successOutcome) * 0, ...
                                successOutcome * 10, ...
                                (isLossState && ~successOutcome) * 0, ...
                                ~successOutcome * 10);
                else
                    reward = fi(isLossState, (successOutcome * 0), (successOutcome * -10), ...
                                (~successOutcome * 0), (~successOutcome * -10));
                end
            end

            if isHighControl
                if obj.outcomeVec(trialIndex) == 1
                    if ~mod(state,2)
                        if action == obj.getCorrectAction(state)
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters (80% of the cases)
                                % if loss state
                                % if action is correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % avoid loss
                                reward = 0;
                            else
                                % if high control
                                % if outcome matters
                                % if loss state
                                % if action is correct
                                % if action leads to incorrect reward
                                % (15% of the cases)
                                % give loss
                                reward = -10;
                            end
                        else
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters
                                % if loss state
                                % if action isn't correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % give loss
                                reward = -10;
                            else
                                % if high control
                                % if outcome matters
                                % if loss state
                                % if action isn't correct
                                % if action leads to incorrect reward
                                % (15% of the cases)
                                % avoid loss
                                reward = 0;
                            end
                        end
                    else
                        if action == obj.getCorrectAction(state)
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters (80% of the cases)
                                % if win state
                                % if action is correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % give reward
                                reward = 10;
                            else
                                % if high control
                                % if outcome matters
                                % if win state
                                % if action is correct
                                % if action leads to incorrect reward
                                % (15% of the cases)
                                % don't give reward
                                reward = 0;
                            end
                        else
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters
                                % if win state
                                % if action isn't correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % don't give reward
                                reward = 0;
                            else
                                % if high control
                                % if outcome matters
                                % if win state
                                % if action isn't correct
                                % if action leads to incorrect reward
                                % (85% of the cases)
                                % give reward
                                reward = 10;
                            end
                        end
                    end
                else
                    if ~mod(state,2)
                        if obj.successVec(trialIndex) == 1
                            % if high control
                            % if outcome doesn't matter (20% of the cases)
                            % if loss state
                            % if action leads to correct reward
                            % (85% of the cases)
                            % avoid loss
                            reward = 0;
                        else
                            % if high control
                            % if outcome matters
                            % if loss state
                            % if action leads to incorrect reward
                            % (15% of the cases)
                            % give loss
                            reward = -10;
                        end
                    else
                        if obj.successVec(trialIndex) == 1
                            % if high control
                            % if outcome doesn't matter (20% of the cases)
                            % if win state
                            % if action leads to correct reward
                            % (85% of the cases)
                            % give loss
                            reward = -10;
                        else
                            % if high control
                            % if outcome doesn't matter
                            % if win state
                            % if action leads to incorrect reward
                            % (85% of the cases)
                            % avoid loss
                            reward = 0;
                        end
                    end
                end
            else
                if obj.outcomeVec(trialIndex) == 0
                    if ~mod(state,2)
                        if action == obj.getCorrectAction(state)
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters (80% of the cases)
                                % if loss state
                                % if action is correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % avoid loss
                                reward = 0;
                            else
                                % if high control
                                % if outcome matters
                                % if loss state
                                % if action is correct
                                % if action leads to incorrect reward
                                % (15% of the cases)
                                % give loss
                                reward = -10;
                            end
                        else
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters
                                % if loss state
                                % if action isn't correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % give loss
                                reward = -10;
                            else
                                % if high control
                                % if outcome matters
                                % if loss state
                                % if action isn't correct
                                % if action leads to incorrect reward
                                % (15% of the cases)
                                % avoid loss
                                reward = 0;
                            end
                        end
                    else
                        if action == obj.getCorrectAction(state)
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters (80% of the cases)
                                % if win state
                                % if action is correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % give reward
                                reward = 10;
                            else
                                % if high control
                                % if outcome matters
                                % if win state
                                % if action is correct
                                % if action leads to incorrect reward
                                % (15% of the cases)
                                % don't give reward
                                reward = 0;
                            end
                        else
                            if obj.successVec(trialIndex) == 1
                                % if high control
                                % if outcome matters
                                % if win state
                                % if action isn't correct
                                % if action leads to correct reward
                                % (85% of the cases)
                                % don't give reward
                                reward = 0;
                            else
                                % if high control
                                % if outcome matters
                                % if win state
                                % if action isn't correct
                                % if action leads to incorrect reward
                                % (85% of the cases)
                                % give reward
                                reward = 10;
                            end
                        end
                    end
                else
                    if ~mod(state,2)
                        if obj.successVec(trialIndex) == 1
                            % if high control
                            % if loss state
                            % if action leads to correct reward
                            % (85% of the cases)
                            % avoid loss
                            reward = 0;
                        else
                            % if high control
                            % if outcome matters
                            % if loss state
                            % if action leads to incorrect reward
                            % (15% of the cases)
                            % give loss
                            reward = -10;
                        end
                    else
                        if obj.successVec(trialIndex) == 1
                            % if high control
                            % if outcome doesn't matter (20% of the cases)
                            % if win state
                            % if action leads to correct reward
                            % (85% of the cases)
                            % give loss
                            reward = -10;
                        else
                            % if high control
                            % if outcome doesn't matter
                            % if win state
                            % if action leads to incorrect reward
                            % (85% of the cases)
                            % avoid loss
                            reward = 0;
                        end
                    end
                end
            end
        end
    end
end