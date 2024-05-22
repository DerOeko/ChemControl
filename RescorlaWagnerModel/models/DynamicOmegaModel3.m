classdef DynamicOmegaModel3 < Model
    % Subclass of Model
    % Omega Model with Qoppa
    properties
        W
        goBias
        V
        Vinit
        omega
        alphaAssoc
        kappaAssoc
        assoc
        assoc0
    end

    methods
        function obj = DynamicOmegaModel3(epsilon, rho, Qinit, goBias, Vinit, alphaAssoc, kappaAssoc, assoc0)
            arguments
                epsilon (1,1) {mustBeNumeric} = config.epsilon
                rho (1,1) {mustBeNumeric} = config.rho
                Qinit (:,:) {mustBeNumeric} = config.Qinit
                goBias (1,1) {mustBeNumeric} = config.goBias
                Vinit (1, :) {mustBeNumeric} = config.Vinit
                alphaAssoc (1, 1) {mustBeNumeric} = config.alphaAssoc
                kappaAssoc (1,1) {mustBeNumeric} = config.kappaAssoc
                assoc0 (1,1) {mustBeNumeric} = config.assoc0
            end
            obj@Model(epsilon, rho, Qinit);

            obj.goBias = goBias;
            obj.alphaAssoc = alphaAssoc;
            obj.kappaAssoc = kappaAssoc;
            obj.Vinit = Vinit;
            obj.V = obj.Vinit;
            obj.assoc0 = assoc0;
            obj.assoc = assoc0;
            obj.omega = 1/(1+exp(-obj.kappaAssoc*(obj.assoc - obj.assoc0)));
            obj.W(:, 1) = (1-obj.omega)*obj.Q(:, 1) + obj.goBias + obj.omega * obj.V(:);
            obj.W(:, 2) = (1-obj.omega)*obj.Q(:,2);
        end

        function obj = calcProbs(obj,state)
            obj.P(state, :) = softmax(obj.W(state, :));
        end

        function obj = updateModel(obj, reward, state, action)

            % Update absolute prediction error
            % Differentiable approximation: 
            %Q_PE = sqrt((reward-obj.Q(state,action)).^2 + eps);
            Q_PE = abs(obj.rho * reward - obj.Q(state, action));

            % Update the Q-values
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * ...
            (obj.rho * reward - obj.Q(state, action));

            % Update V-values
            obj.V(state) = obj.V(state) + obj.epsilon * ...
            (obj.rho * reward - obj.V(state));
            
            obj.assoc = obj.assoc + obj.alphaAssoc * obj.epsilon * (Q_PE - obj.assoc);
            disp(obj.assoc);
            %obj.omega = 1/(1+exp(-obj.kappaAssoc*(obj.assoc - obj.assoc0)));
            obj.omega = min(obj.kappaAssoc*(obj.assoc-obj.assoc0), 1);

            if action == 1
                obj.W(state, 1) = (1-obj.omega)*obj.Q(state, 1) + obj.goBias + obj.omega * obj.V(state);
            else
                obj.W(state, 2) = (1-obj.omega)*obj.Q(state, 2);
            end
        end

        function obj = resetQ(obj)
            obj.Q = obj.Qinit;
            obj.V = obj.Vinit;
            obj.W(:, 1) = (1-obj.omega)*obj.Q(:, 1) + obj.goBias  + obj.omega * obj.V(:);
            obj.W(:, 2) = (1-obj.omega)*obj.Q(:,2);
        end

        function omega = getOmega(obj)
            omega = obj.omega;
        end
    end
end