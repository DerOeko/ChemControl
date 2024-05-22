classdef FixedOmegaModel < DynamicMotivationalBiasModel
    % Subclass of DynamicMotivationalBiasModel with fixed omega
    properties
        o
    end

    methods
        function obj = FixedOmegaModel(cfg)
            arguments
                cfg Config
            end
            obj@DynamicMotivationalBiasModel(cfg);

            obj.o = cfg.o;
        end

        function obj = calcWeights(obj, s)
            obj.W(s, 1) = (1-obj.o) * obj.Q(s, 1) + obj.gB + obj.o * obj.V(s);
            obj.W(s, 2) = (1-obj.o) * obj.Q(s, 2);
        end

        function o = getOmega(obj)
            o = obj.o;
        end
    end
end

