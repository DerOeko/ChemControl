classdef DynamicMotivationalBiasModel < GoBiasModel
    % Subclass of Model with GoBiasModel
    properties
        Vinit
        V
        pi
    end

    methods
        function obj = DynamicMotivationalBiasModel(epsilon, rho, beta, Qinit, goBias, Vinit, pi)
            % Define default values using arguments block
            arguments
                epsilon
                rho
                beta
                Qinit = zeros(4,2)
                goBias = 0.5
                Vinit = [0.5 -0.5 0.5 -0.5] % Default value for V
                pi = 0.5 % Default value for pi
            end
            obj@GoBiasModel(epsilon, rho, beta, Qinit, goBias);

            % Check that V has length matching Q's number of rows
            if length(Vinit) ~= size(obj.Q, 1)
                error('Length of Vinit must match the number of states in Q.');
            end
        
            obj.Vinit = Vinit;
            obj.V = obj.Vinit;
            obj.pi = pi;
            obj.W(:, 1) = obj.W(:,1) + obj.Vinit(:)*obj.pi;
        end

        function obj = updateModel(obj, reward, state, action)
          % Update the Q-values
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
            obj.V(state) = obj.V(state) + obj.epsilon*(obj.rho* - obj.V(state));
            % Update the W matrix
            % For the "Go" action, update the value with the Q-value +
            % goBias + fixed motivational bias.
            obj.W(state, 1) = obj.Q(state, 1) + obj.goBias + obj.V(state)*obj.pi;
            % For the "NoGo" action, update the value with the Q-value
            obj.W(state, 2) = obj.Q(state, 2);
        end

        function obj = resetQ(obj)
            obj.Q = obj.Qinit;
            obj.W(:, 1) = obj.Q(:,1) + obj.goBias + obj.V(:)*obj.pi;
            obj.W(:,2) = obj.Q(:,2);
            obj.V = obj.Vinit;
        end
    end
end

