classdef Model2
    properties
        Q
        P
        epsilon
        beta
        rho
        Qinit
    end
    methods
        function obj = Model2(epsilon, rho, beta, Qinit)
            obj.epsilon = epsilon;
            obj.rho = rho;
            obj.beta = beta;
            obj.Qinit = Qinit;
            obj.Q = obj.Qinit;
            obj.P = zeros(4,2)+0.5;
        end

        function obj = calcProbs(obj, state)
            obj.P(state, :) = betaSoftmax(obj.Q(state, :), obj.beta);
        end

        function action = returnAction(obj, state)
            action = randsample([1,2], 1, true, obj.P(state, :));
        end

        function obj = updateModel(obj, reward, state, action)
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * ( obj.rho * reward - obj.Q(state,action));
        end

        function obj = resetQ(obj)
            obj.Q = obj.Qinit;
        end

        function obj = resetP(obj)
            obj.P = zeros(4,2)+0.5;
        end

        function goProb = returnGoProb(obj, state)
            goProb = obj.P(state, 1);
        end
    end
end
