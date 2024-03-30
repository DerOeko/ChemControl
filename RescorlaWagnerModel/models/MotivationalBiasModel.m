classdef MotivationalBiasModel < GoBiasModel
    % Subclass of Model with Go Bias
    properties
        V
        pi
    end

    methods
        function obj = MotivationalBiasModel(epsilon, rho, beta, goBias, V, pi)
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

            % Check that goBias is a scalar
            if ~isscalar(goBias)
                error('goBias must be a scalar.');
            end
            
            % Check that V has length matching Q's number of rows
            if length(V) ~= size(obj.Q, 1)
                error('Length of V must match the number of states in Q.');
            end

            obj.V = V;
            obj.pi = pi;
            obj.W = obj.Q;
            obj.W(:, 1) = obj.W(:,1) + obj.V(:)*obj.pi;
        end

        function obj = updateModel(obj, reward, state, action)
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
            if action == 1
                obj.W(state, action) = obj.Q(state,action) + obj.goBias + obj.V(state)*obj.pi;
            else
                obj.W(state,action) =obj.Q(state,action);
            end
        end

        function obj = resetQ(obj)
            obj.Q = ones(4,2);
            obj.W(:, 1) = obj.W(:,1) + obj.goBias + obj.V(:)*obj.pi;
            obj.W(:,2) = obj.Q(:,2);
        end
    end
end

