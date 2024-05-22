classdef TrialEnvironment
    properties
        % Probability of the correct action leading to a good outcome
        rP = 0.85;
        
        % Probability of the action mattering
        cP = 0.8;
        
        % Assuming env.cond = {'GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss'}
        % and actions are represented as [1, 2] for ['Go', 'NoGo']
        cond = {1, 2, 3, 4};
        
        % Array of states
        sArr;
        
        % Number of trials in a block
        T;
        
        % Vector indicating whether the outcome matters for each trial
        cVec;
        
        % Vector indicating whether the action leads to the correct reward/penalty for each trial
        rVec;
    end

    methods
        function obj = TrialEnvironment(rP, cP, sArr, T)
            % Constructor
            % Set the probability of the correct action leading to a good outcome
            obj.rP = rP;
            
            % Set the probability of the action mattering
            obj.cP = cP;
            
            % Set the array of states
            obj.sArr = sArr;

            % Calculate the number of trials where the outcome matters
            O = round(cP * T);
            
            % Generate a vector with 1s for trials where the outcome matters,
            % and 0s for trials where the outcome doesn't matter
            obj.cVec = [ones(1, O) zeros(1, T - O, 1)];

            % Shuffle the cVec
            obj.cVec = obj.cVec(randperm(T));
            
            % Calculate the number of trials where the action leads to the correct reward/penalty
            S = round(rP * T);
            
            % Generate a vector with 1s for trials where the action leads to the correct reward/penalty,
            % and 0s for trials where the action doesn't lead to the correct reward/penalty
            obj.rVec = [ones(1, S) zeros(1, T - S, 1)];

            % Shuffle the rVec
            obj.rVec = obj.rVec(randperm(T));
        end
        
        function [state, correctAction, obj] = presentTrial(obj)
            % Present a trial
            % Get the current state from the sArr
            state = obj.sArr(1);
            
            % Remove the current state from the sArr
            obj.sArr(1) = [];
            
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
        
        function reward = getReward(obj, trialIndex, state, action, isHC)
            % Calculate the reward for a given trial
            
            % Get whether the outcome matters for the current trial
            outcomeMatters = obj.cVec(trialIndex);
            
            % Get whether the action leads to the correct reward/penalty for the current trial
            successOutcome = obj.rVec(trialIndex);
            
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

            if isHC
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