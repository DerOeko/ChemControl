classdef DynamicOmegaModel1 < FixedOmegaModel
    % Subclass of Model
    % Omega Model with Qoppa
    properties
        o0
        alpha
        kap
        slope
        alphaQ
        alphaV
    end

    methods
        function obj = DynamicOmegaModel1(cfg)
            arguments
                cfg Config
            end
            obj@FixedOmegaModel(cfg);


            obj.o0 = cfg.o0;
            obj.o = obj.o0;
            obj.alpha = cfg.alpha;
            obj.kap = cfg.kap;
            obj.slope = cfg.slope;

            obj.alphaQ = 0;
            obj.alphaV = 0;

        end

        function obj = updateModel(obj, r, s, a)

            % Update absolute prediction errors
            Q_PE = abs(obj.rho * r - obj.Q(s, a));
            V_PE = abs(obj.rho * r - obj.V(s));

            % Update the Q-values
            obj.Q(s, a) = obj.Q(s, a) + obj.ep * (obj.rho * r - obj.Q(s, a));

            % Update V-values
            obj.V(s) = obj.V(s) + obj.ep * (obj.rho * r - obj.V(s));
            
            % Update "Q" and "V" trust
            obj.alphaQ = obj.alphaQ + obj.alpha * obj.ep * (Q_PE - obj.alphaQ);
            obj.alphaV = obj.alphaV + obj.alpha * obj.ep * (V_PE - obj.alphaV);

            % Calculate relative ratio between uncertainty of Q and V. If Q
            % is more uncertain, be bigger.
            ratio = obj.alphaQ/(obj.alphaQ + obj.alphaV);
            
            % Update o and soft clip to [0, 1]
            obj.o = 1/(1 + exp(-obj.slope * ((obj.o + obj.kap * (ratio - 0.5))-0.5)));
        end
    end
end