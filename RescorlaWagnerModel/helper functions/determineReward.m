function reward = determineReward(outcome, success, state, correctAction, isHighControl)
    % Calculate the reward for a given trial
            
    % Get whether the outcome matters for the current trial
    outcomeMatters = outcome;
    
    % Get whether the action leads to the correct reward/penalty for the current trial
    successOutcome = success;
    
    % Check if the current state is a loss state or not
    isLossState = ~mod(state, 2);

    fi = @(varargin) varargin{length(varargin) - varargin{1}};

    if isHighControl
        % If it's a high control block
        if outcomeMatters
            % If the outcome matters
            if correctAction
                % If the agent's action is correct
                if isLossState
                    % If it's a loss state
                    reward = fi(successOutcome, 0, -10);
                else
                    % If it's a win state
                    reward = fi(successOutcome, 10, 0);
                end
            else
                % If the agent's action is incorrect
                if isLossState
                    % If it's a loss state
                    reward = fi(~successOutcome, 0, -10);
                else
                    % If it's a win state
                    reward = fi(~successOutcome, 10, 0);
                end
            end
        else
            % If the outcome doesn't matter
            if isLossState
                % If it's a loss state
                reward = fi(successOutcome, 0, -10);
            else
                % If it's a win state
                reward = fi(successOutcome, 10, 0);
            end
        end
    else
        % If it's not a high control block
        if ~outcomeMatters
            % If the outcome doesn't matter
            if correctAction
                % If the agent's action is correct
                if isLossState
                    % If it's a loss state
                    reward = fi(successOutcome, 0, -10);
                else
                    % If it's a win state
                    reward = fi(successOutcome, 10, 0);
                end
            else
                % If the agent's action is incorrect
                if isLossState
                    % If it's a loss state
                    reward = fi(~successOutcome, 0, -10);
                else
                    % If it's a win state
                    reward = fi(~successOutcome, 10, 0);
                end
            end
        else
            % If the outcome matters
            if isLossState
                % If it's a loss state
                reward = fi(successOutcome, 0, -10);
            else
                % If it's a win state
                reward = fi(successOutcome, 10, 0);
            end
        end
    end
end