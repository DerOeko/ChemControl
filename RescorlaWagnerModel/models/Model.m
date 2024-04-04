classdef Model
    % Super class of basic rescorla wagner model
    properties
        Q
        Qinit
        P
        epsilon
        beta
        rho
    end

    methods
        function obj = Model(epsilon, rho, beta, Qinit)
            arguments
                % Initialize learning rate (epsilon) for Q-value updates
                epsilon = 0.31
                % Feedback sensitivity (rho) scales the difference between received and estimated reward
                rho = 5
                % Initialize softmax temperature (beta) for the softmax action selection
                beta = 3

                Qinit = zeros(4,2);
            end
            obj.epsilon = epsilon;
            obj.rho = rho;
            obj.beta = beta;
            obj.Qinit = Qinit;
            obj.Q = obj.Qinit;
            obj.P = zeros(4,2) + 0.5;
        end

        function obj = calcProbs(obj,state)
            obj.P(state, :) = betaSoftmax(obj.Q(state, :), obj.beta);
        end

        function action = returnAction(obj, state)
            action = randsample([1,2], 1, true, obj.P(state, :));
        end

        function obj = updateModel(obj, reward, state, action)
            obj.Q(state, action) = obj.Q(state, action) + obj.epsilon * (obj.rho * reward - obj.Q(state, action));
        end

        function obj = resetQ(obj)
            obj.Q = obj.Qinit;
        end

        function obj = resetP(obj)
            obj.P = zeros(4,2) + 0.5;
        end

        function goProb = returnGoProb(obj, state)
            goProb = obj.P(state, 1);
        end
    end
end

