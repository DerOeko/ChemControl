classdef Model
    % Super class of basic rescorla wagner model
    properties
        Q
        Q0
        P
        ep
        rho
        W
    end

    methods
        function obj = Model(cfg)
            arguments
                cfg Config
            end
            
            obj.ep = cfg.ep;
            obj.rho = cfg.rho;
            obj.Q0 = cfg.Q0;
            obj.Q = obj.Q0;
            obj.P = zeros(4,2) + 0.5;
        end

        function obj = calcProbs(obj, s)
            obj.P(s, 1) = 1./(1+exp(-(obj.W(s, 1)-obj.W(s, 2))));
            obj.P(s, 2) = 1- obj.P(s,1); 
        end

        function a = returnAction(obj, s)
            a = randsample([1,2], 1, true, obj.P(s, :));
        end
        
        function obj = calcWeights(obj, s)
            obj.W(s, 1) = obj.Q(s, 1);
            obj.W(s, 2) = obj.Q(s, 2);
        end

        function obj = updateModel(obj, r, s, a)
            obj.Q(s, a) = obj.Q(s, a) + obj.ep * ((obj.rho * r) - obj.Q(s, a));
        end

        function obj = resetModel(obj)
            obj.Q = obj.Q0;
            obj.P = zeros(4,2) + 0.5;

        end

        function goProb = returnGoProb(obj, s)
            goProb = obj.P(s, 1);
        end
    end
end
