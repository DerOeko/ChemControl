classdef GoBiasModel < Model
    % Subclass of Model with Go Bias
    properties
        W
        goBias
    end

    methods
        function obj = GoBiasModel(epsilon, rho, beta, goBias)
            arguments
                epsilon
                rho
                beta
                goBias = 0.5
            end
            obj@Model(epsilon, rho, beta);

            % Check that goBias is a scalar
            if ~isscalar(goBias)
                error('goBias must be a scalar.');
            end

            obj.goBias = goBias;
            obj.W = obj.Q;
            obj.W(:, 1) = obj.W(:,1) + obj.goBias;
            disp(obj.W)
        end

        function obj = calcProbs(obj,state)
            obj.P(state, :) = betaSoftmax(obj.W(state, :), obj.beta);
        end

        function obj = updateModel(obj, reward, state, action)
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
            if action == 1
                obj.W(state, action) = obj.Q(state,action) + obj.goBias;
            else
                obj.W(state,action) =obj.Q(state,action);
            end
        end

        function obj = resetQ(obj)
            obj.Q = ones(4,2);
            obj.W(:, 1) = obj.W(:,1) + obj.goBias;
            obj.W(:,2) = obj.Q(:,2);        
        end
    end
end

