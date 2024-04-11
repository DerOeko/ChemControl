classdef TrialEnvironment
    properties
        % Probability of the correct action leading to a good outcome
        successProb = 0.85;
        
        % Probability of the action mattering
        outcomeProb = 0.8;
        
        % Assuming env.conditions = {'GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss'}
        % and actions are represented as [1, 2] for ['Go', 'NoGo']
        conditions = {1, 2, 3, 4};
        
        % Array of states
        stateArray;
        
        % Number of trials in a block
        numTrialsInBlock;
        
        % Vector indicating whether the outcome matters for each trial
        outcomeVec;
        
        % Vector indicating whether the action leads to the correct reward/penalty for each trial
        successVec;
    end

    methods
        function obj = TrialEnvironment(successProb, outcomeProb, stateArray, numTrialsInBlock)
            % Constructor
            % Set the probability of the correct action leading to a good outcome
            obj.successProb = successProb;
            
            % Set the probability of the action mattering
            obj.outcomeProb = outcomeProb;
            
            % Set the array of states
            obj.stateArray = stateArray;

            % Calculate the number of trials where the outcome matters
            nOutcomes = round(outcomeProb * numTrialsInBlock);
            
            % Generate a vector with 1s for trials where the outcome matters,
            % and 0s for trials where the outcome doesn't matter
            obj.outcomeVec = [ones(1, nOutcomes) zeros(1, numTrialsInBlock - nOutcomes, 1)];

            % Shuffle the outcomeVec
            obj.outcomeVec = obj.outcomeVec(randperm(numTrialsInBlock));
            
            % Calculate the number of trials where the action leads to the correct reward/penalty
            nSuccesses = round(successProb * numTrialsInBlock);
            
            % Generate a vector with 1s for trials where the action leads to the correct reward/penalty,
            % and 0s for trials where the action doesn't lead to the correct reward/penalty
            obj.successVec = [ones(1, nSuccesses) zeros(1, numTrialsInBlock - nSuccesses, 1)];

            % Shuffle the successVec
            obj.successVec = obj.successVec(randperm(numTrialsInBlock));
        end
        
        function [state, correctAction, obj] = presentTrial(obj)
            % Present a trial
            % Get the current state from the stateArray
            state = obj.stateArray(1);
            
            % Remove the current state from the stateArray
            obj.stateArray(1) = [];
            
            % Get the correct action for the current state
            correctAction = obj.getCorrectAction(state);
        end
        
        function correctAction = getCorrectAction(obj, state)
            % Get the correct action for a given state
            if state <= 2
                correctAction = 1; % 'Go'
            else
                correctAction = 2; % 'NoGo'
            end
        end
        
        function reward = getReward(obj, trialIndex, state, action, isHighControl)
            % Calculate the reward for a given trial
            
            % Get whether the outcome matters for the current trial
            outcomeMatters = obj.outcomeVec(trialIndex);
            
            % Get whether the action leads to the correct reward/penalty for the current trial
            successOutcome = obj.successVec(trialIndex);
            
            % Check if the current state is a loss state or not
            isLossState = ~mod(state, 2);
            
            % Check if the agent's action is correct or not
            correctAction = action == obj.getCorrectAction(state);
            
            % Define the fi function to model a ternary operator from Python/Java
            % The fi function takes a condition as the first argument,
            % and then pairs of values as the remaining arguments.
            % It returns the value corresponding to the first condition that evaluates to true.
            %
            % For example, fi(condition, value1, value2) will return:
            %   - value1 if condition is true (non-zero value)
            %   - value2 if condition is false (zero)
            %
            % The fi function works as follows:
            %   1. The first argument (varargin{1}) is treated as a condition.
            %   2. If the condition is true (non-zero value), it returns the last argument (varargin{length(varargin)}).
            %   3. If the condition is false (zero), it returns the second-to-last argument (varargin{length(varargin)-1}).
            fi = @(varargin) varargin{length(varargin) - varargin{1}};

            if isHighControl
                % If it's a high control block
                if outcomeMatters
                    % If the outcome matters
                    if correctAction
                        % If the agent's action is correct
                        if isLossState
                            % If it's a loss state
                            reward = fi(successOutcome, 0, -1);
                        else
                            % If it's a win state
                            reward = fi(successOutcome, 1, 0);
                        end
                    else
                        % If the agent's action is incorrect
                        if isLossState
                            % If it's a loss state
                            reward = fi(~successOutcome, 0, -1);
                        else
                            % If it's a win state
                            reward = fi(~successOutcome, 1, 0);
                        end
                    end
                else
                    % If the outcome doesn't matter
                    if isLossState
                        % If it's a loss state
                        reward = fi(successOutcome, 0, -1);
                    else
                        % If it's a win state
                        reward = fi(successOutcome, 1, 0);
                    end
                end
            else
                % If it's not a high control block
                if ~outcomeMatters
                    % In 20% of the cases, action matters
                    if correctAction
                        % If the agent's action is correct
                        if isLossState
                            % If it's a loss state
                            reward = fi(successOutcome, 0, -1);
                        else
                            % If it's a win state
                            reward = fi(successOutcome, 1, 0);
                        end
                    else
                        % If the agent's action is incorrect
                        if isLossState
                            % If it's a loss state
                            reward = fi(~successOutcome, 0, -1);
                        else
                            % If it's a win state
                            reward = fi(~successOutcome, 1, 0);
                        end
                    end
                else
                    % If the outcome matters
                    if isLossState
                        % If it's a loss state
                        reward = fi(successOutcome, 0, -1);
                    else
                        % If it's a win state
                        reward = fi(successOutcome, 1, 0);
                    end
                end
            end
        end
    end
end