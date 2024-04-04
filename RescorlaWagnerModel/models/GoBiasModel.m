classdef GoBiasModel < Model
    % Subclass of Model with Go Bias
    properties
        W
        goBias
    end

    methods
        function obj = GoBiasModel(epsilon, rho, beta, Qinit, goBias)
            arguments
                epsilon
                rho
                beta
                Qinit = zeros(4,2)
                goBias = 0.3
            end
            obj@Model(epsilon, rho, beta, Qinit);

            % Check that goBias is a scalar
            if ~isscalar(goBias)
                error('goBias must be a scalar.');
            end

            obj.goBias = goBias;
            obj.W(:, 1) = obj.Q(:,1) + obj.goBias;
            obj.W(:, 2) = obj.Q(:,2);
        end

        function obj = calcProbs(obj,state)
            obj.P(state, :) = betaSoftmax(obj.W(state, :), obj.beta);
        end

        function obj = updateModel(obj, reward, state, action)
            % Update the Q-values
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
        
            % Update the W matrix
            % For the "Go" action, update the value with the Q-value + goBias
            obj.W(state, 1) = obj.Q(state, 1) + obj.goBias;
            % For the "NoGo" action, update the value with the Q-value
            obj.W(state, 2) = obj.Q(state, 2);
        end

        function obj = resetQ(obj)
            obj.Q = obj.Qinit;
            obj.W(:, 1) = obj.Q(:,1) + obj.goBias;
            obj.W(:,2) = obj.Q(:,2);
        end
    end
end

