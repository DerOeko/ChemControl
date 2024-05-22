classdef FixedMotivationalBiasModel < GoBiasModel
    % Subclass of Model with Go Bias
    properties
        V
        pi
    end

    methods
        function obj = FixedMotivationalBiasModel(cfg)
            arguments
                cfg Config
            end
            obj@GoBiasModel(cfg);

            obj.V = cfg.V;
            obj.pi = cfg.pi;
        end

        function obj = calcWeights(obj, s)
            obj.W(s, 1) = obj.Q(s, 1) + obj.gB + obj.pi * obj.V(s);
            obj.W(s, 2) = obj.Q(s, 2);
        end
    end
end

