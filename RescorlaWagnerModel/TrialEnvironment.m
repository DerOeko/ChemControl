classdef TrialEnvironment
    properties
        successProb = 0.85;
        outcomeProb = 0.8;
        % Assuming env.conditions = {'GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss'}
        % and actions are represented as [1, 2] for ['Go', 'NoGo']
        conditions = {1,2,3,4};
        stateArray;
    end

    methods
        function obj = TrialEnvironment(successProb, outcomeProb, stateArray)
            % How likely is it that the correct action leads to a good
            % outcome?
            obj.successProb = successProb;
            % How likely is it that the action matters?
            obj.outcomeProb = outcomeProb;

            obj.stateArray = stateArray;
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

        function reward = getReward(obj, state, action, isHighControl)
            % If it is a high control block, then in 80% or w/e of the
            % cases, your actions matter.
            if isHighControl
                % This is true in obj.outcomeProb% of the cases.
                if rand < obj.outcomeProb
                    % This looks at whether it is loss or win.
                    if ~mod(state, 2)
                        % This is loss
                        if action == obj.getCorrectAction(state)
                            % In 85% percent of the cases, the bool is 0.
                            % Thus, you get 0 reward (or avoid loss).
                            reward = -10 * (rand >= obj.successProb);
                        else
                            reward = -10 * (rand < obj.successProb);
                        end
                    % Win
                    else
                        if action == obj.getCorrectAction(state)
                            % In 85% percent of the cases, the bool is 1.
                            % Then, you get 10 reward. Otherwise you get 0.
                            reward = 10 * (rand < obj.successProb);
                        else
                            reward = 10 * (rand >= obj.successProb);
                        end
                    end
                % Even in a high control block, there are trials where the
                % outcomes do not matter.
                else
                    if ~mod(state, 2)
                        reward = -10 * (rand >= obj.successProb);
                    else
                        reward = 10 * (rand < obj.successProb);
                    end
                end
            else
                % Here, it is inversed. In 20% of the cases, this is true.
                if rand < 1 - obj.outcomeProb
                    % Now, actions matter, even though ~isHighControl.
                    if ~mod(state, 2)
                        if action == obj.getCorrectAction(state)
                            % In 85% percent of the cases, the bool is 0.
                            % Thus, you get 0 reward (or avoid loss).
                            reward = -10 * (rand >= obj.successProb);
                        else
                            % Check this line as well
                            reward = -10 * (rand < obj.successProb);
                        end
                    else
                        if action == obj.getCorrectAction(state)
                            reward = 10 * (rand < obj.successProb);
                        else
                            reward = 10 * (rand >= obj.successProb);
                        end
                    end
                else
                    if ~mod(state, 2)
                        reward = -10 * (rand >= obj.successProb);
                    else
                        reward = 10 * (rand < obj.successProb);
                    end
                end
            end
        end
    end
end