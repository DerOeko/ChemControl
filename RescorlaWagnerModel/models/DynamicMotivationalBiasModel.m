classdef DynamicMotivationalBiasModel < FixedMotivationalBiasModel
    % Subclass of Model with GoBiasModel
    properties
        V0
    end

    methods
        function obj = DynamicMotivationalBiasModel(cfg)
            arguments
                cfg Config
            end
            obj@FixedMotivationalBiasModel(cfg);
            obj.V0 = cfg.V0;
            obj.V = obj.V0;
        end

        function obj = updateModel(obj, r, s, a)
            obj.Q(s, a) = obj.Q(s, a) + obj.ep * ((obj.rho * r) - obj.Q(s, a));
            obj.V(s) = obj.V(s) + obj.ep*(obj.rho * r - obj.V(s));
        end

        function obj = resetModel(obj)
            obj.Q = obj.Q0;
            obj.P = zeros(4,2) + 0.5;
            obj.V = obj.V0;
        end
    end
end

