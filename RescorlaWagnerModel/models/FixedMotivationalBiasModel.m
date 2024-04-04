classdef FixedMotivationalBiasModel < GoBiasModel
    % Subclass of Model with Go Bias
    properties
        V
        pi
    end

    methods
        function obj = FixedMotivationalBiasModel(epsilon, rho, beta, goBias, V, pi)
            % Define default values using arguments block
            arguments
                epsilon
                rho
                beta
                goBias = 0.5
                V = [0.5 -0.5 0.5 -0.5] % Default value for V
                pi = 0.5 % Default value for pi
            end
            obj@GoBiasModel(epsilon, rho, beta, goBias);

            % Check that V has length matching Q's number of rows
            if length(V) ~= size(obj.Q, 1)
                error('Length of V must match the number of states in Q.');
            end

            obj.V = V;
            obj.pi = pi;
            obj.W(:, 1) = obj.W(:,1) + obj.V(:)*obj.pi;
            disp(obj.W)
        end

        function obj = updateModel(obj, reward, state, action)
          % Update the Q-values
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
        
            % Update the W matrix
            % For the "Go" action, update the value with the Q-value +
            % goBias + fixed motivational bias.
            obj.W(state, 1) = obj.Q(state, 1) + obj.goBias + obj.V(state)*obj.pi;
            % For the "NoGo" action, update the value with the Q-value
            obj.W(state, 2) = obj.Q(state, 2);
        end

        function obj = resetQ(obj)
            obj.Q = zeros(4,2) + 0.5;
            obj.W(:, 1) = obj.Q(:,1) + obj.goBias + obj.V(:)*obj.pi;
            obj.W(:,2) = obj.Q(:,2);
        end
    end
end

