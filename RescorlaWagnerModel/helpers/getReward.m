function o = getReward(t, cVec, rVec, state, action, isHC)
    % Calculate the o for a given trial
    
    % Get whether the outcome matters for the current trial
    outcomeMatters = cVec(t);
    
    % Get whether the action leads to the correct o/penalty for the current trial
    successOutcome = rVec(t);
    
    % Check if the current state is a loss state or not
    isLossState = ~mod(state, 2);
    
    % Check if the agent's action is correct or not
    correctAction = action ==  getCorrectAction(state);
    
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
                    o = fi(successOutcome, 0, -1);
                else
                    % If it's a win state
                    o = fi(successOutcome, 1, 0);
                end
            else
                % If the agent's action is incorrect
                if isLossState
                    % If it's a loss state
                    o = fi(~successOutcome, 0, -1);
                else
                    % If it's a win state
                    o = fi(~successOutcome, 1, 0);
                end
            end
        else
            % If the outcome doesn't matter
            if isLossState
                % If it's a loss state
                o = fi(successOutcome, 0, -1);
            else
                % If it's a win state
                o = fi(successOutcome, 1, 0);
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
                    o = fi(successOutcome, 0, -1);
                else
                    % If it's a win state
                    o = fi(successOutcome, 1, 0);
                end
            else
                % If the agent's action is incorrect
                if isLossState
                    % If it's a loss state
                    o = fi(~successOutcome, 0, -1);
                else
                    % If it's a win state
                    o = fi(~successOutcome, 1, 0);
                end
            end
        else
            % If the outcome matters
            if isLossState
                % If it's a loss state
                o = fi(successOutcome, 0, -1);
            else
                % If it's a win state
                o = fi(successOutcome, 1, 0);
            end
        end
    end
end