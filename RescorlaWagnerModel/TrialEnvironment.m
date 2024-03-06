classdef TrialEnvironment
    properties
        successProbability = 0.85;
        conditions = {'GoToWin', 'GoToAvoidLoss', 'NoGoToWin', 'NoGoToAvoidLoss'}
    end

    methods
        function obj = TrialEnvironment()
        end
        
        function [state, correctAction] = presentTrial(obj)
            state = randsample(obj.conditions, 1);
            
            correctAction = obj.getCorrectAction(state);
        end

        function correctAction = getCorrectAction(obj, state)
            if any(strcmp(state, {'GoToWin', 'GoToAvoidLoss'}))
                correctAction = 'Go';
            else
                correctAction = 'NoGo';
            end
        end

        function reward = getReward(obj, state, action)
            if endsWith(state, 'Loss') 
                if strcmp(action, obj.getCorrectAction(state))
                    if rand < obj.successProbability
                        reward = 0;
                    else
                        reward = -10;
                    end
                else
                    reward = -10;
                end
            else
                if strcmp(action, obj.getCorrectAction(state))
                    if rand < obj.successProbability
                        reward = 10;
                    else
                        reward = 0;
                    end
                else
                    reward = 0;
                end
            end
        end
    end
end