classdef Model
    properties
        Q
        P
        epsilon
        beta
        rho
    end

    methods
        function obj = Model(epsilon, rho, beta)
            obj.epsilon =epsilon;
            obj.rho = rho;
            obj.beta = beta;
            obj.Q = zeros(4,2) +1;
            obj.P = zeros(4,2) +0.5;
        end

        function obj = calcProbs(obj,state)
            obj.P(state, :) = beta_softmax(obj.Q(state, :), obj.beta);
        end

        function action = returnAction(obj, state)
            action = randsample([1,2], 1, true, obj.P(state, :));
        end

        function obj = updateModel(obj, reward, state, action)
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
        end

        function obj = resetQ(obj)
            obj.Q = zeros(4,2) +1;
        end

        function obj = resetP(obj)
            obj.P = zeros(4,2) + 0.5;
        end

        function goProb = returnGoProb(obj, state)
            goProb = obj.P(state, 1);
        end
    end
end

