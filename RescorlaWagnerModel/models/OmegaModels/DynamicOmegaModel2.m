classdef DynamicOmegaModel2 < FixedOmegaModel
    % Subclass of Model
    % Om Model with Om
    properties
        o0
        Om
        Om0
        aO
        bO
        tO
    end

    methods
        function obj = DynamicOmegaModel2(cfg)
            arguments
                cfg Config
            end
            obj@FixedOmegaModel(cfg);

            obj.o0 = cfg.o0;
            obj.o = obj.o0;
            obj.aO = cfg.aO;
            obj.bO = cfg.bO;
            obj.tO = cfg.tO;
            obj.Om0 = cfg.Om0;
            obj.Om = obj.Om0;
        end

        function obj = updateModel(obj, r, s, a)

            % Calculcate signed prediction errors
            Q_PE = obj.rho * r - obj.Q(s, a);
            V_PE = obj.rho * r - obj.V(s);
            
            % Update the Q-values
            obj.Q(s, a) = obj.Q(s, a) + obj.ep * (obj.rho * r - obj.Q(s, a));

            % Update V-values
            obj.V(s) = obj.V(s) + obj.ep * (obj.rho * r - obj.V(s));
            
            obj.Om = obj.Om + obj.aO*(Q_PE - V_PE - obj.Om);

            obj.o = 1/(1+exp(-obj.bO*(obj.Om-obj.tO)));
        end
    end
end

