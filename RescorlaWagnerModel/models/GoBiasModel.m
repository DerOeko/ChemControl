classdef GoBiasModel < Model
    % Subclass of Model with Go Bias
    properties
        gB
    end

    methods
        function obj = GoBiasModel(cfg)
            arguments
                cfg Config
            end
            obj@Model(cfg);

            obj.gB = cfg.gB;
        end

        function obj = calcWeights(obj, s)
            obj.W(s, 1) = obj.Q(s, 1) + obj.gB;
            obj.W(s, 2) = obj.Q(s, 2);
        end

    end
end
