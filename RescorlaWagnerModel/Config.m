classdef Config
    properties
        % Standard Parameters
        T % Number of Trials in Block
        B % Number of Blocks
        rP % Probability of receiving a reward.
        cP % Probability that the action matters in determining whether to reward or not.
        HC_B % Number of high control blocks
        LC_B
        R % Number of Runs
        cArr % Array of 1's or 0's indicating whether the current block is high or low control.
        
        % Basic Model Parameters
        ep % epsilon
        rho % rho
        Q0 % Qinit
        
        % GoBiasModel Parameters
        gB % goBias
        
        % FixedMotivationalBiasModel Parameters
        V
        pi
        
        % DynamicMotivationalBiasModel Parameters
        V0 % Vinit
        
        % FixedOmegaModel Parameters
        o % omega
        
        % DynamicOmegaModel with Qoppa Parameters
        o0 % omegaInit
        alpha % scaling factor for 
        kap % kappa
        slope
        
        % DynamicOmegaModel with Omega Parameters
        Om0 % OmegaInit
        aO % alphaOmega
        bO % betaOmega
        tO % thresOmega
        
        % DynamicOmegaModel with Associability Parameters
        assoc0
        kA % kappaAssoc
        aA % alphaAssoc
        rA % rhoAssoc
        
        % Plotting Parameters
        lS % lineStyles
        model_names
        omegaModelNames
        models
        numOmegas
    end
    
    methods
        function obj = Config()
            % Initialize standard Parameters
            obj.T = 160;
            obj.B = 16;
            obj.rP = 0.85;
            obj.cP = 0.8;
            obj.HC_B = obj.B / 2;
            obj.LC_B = obj.B - obj.HC_B;
            obj.R = 5;
            
            % Initialize controllabilityArray
            obj.cArr = zeros(obj.B, 1);
            for i = 1:length(obj.cArr)
                if mod(i, 2)
                    obj.cArr(i) = 1;
                else
                    obj.cArr(i) = 0;
                end
            end
            
            % Initialize Model Parameters
            obj.ep = 0.3;
            obj.rho = 4;
            obj.Q0 = reshape([0.5 -0.5 0.5 -0.5, 0.5 -0.5 0.5 -0.5] * obj.rho, [4,2]);
            
            % Initialize GoBiasModel Parameters
            obj.gB = 0.3108;
            
            % Initialize FixedMotivationalBiasModel Parameters
            obj.V = [0.5 -0.5 0.5 -0.5] * obj.rho;
            obj.pi = 0.3;
            
            % Initialize DynamicMotivationalBiasModel Parameters
            obj.V0 = obj.V;
            
            % Initialize FixedOmegaModel Parameters
            obj.o = 0.3;
            
            % Initialize DynamicOmegaModel with Qoppa Parameters
            obj.o0 = 0.3;
            obj.alpha = 0.2;
            obj.kap = 3;
            obj.slope = 3;
            
            % Initialize DynamicOmegaModel with Omega Parameters
            obj.Om0 = 0;
            obj.aO = 0.2;
            obj.bO = 5;
            obj.tO = 0.6;
            
            % Initialize DynamicOmegaModel with Associability Parameters
            obj.assoc0 = 0;
            obj.kA = 1;
            obj.aA = 0.1;
            obj.rA = 1;
            
            % Initialize Plotting Parameters
            obj.lS = {'-', '--', ':', '-.'};

            obj.model_names = {'Generic Model', 'Go Bias Model', 'Fixed Motivational Bias Model', 'Dynamic Motivational Bias Model', "Fixed Omega Model", "Dynamic Omega Model with Qoppa", "Dynamic Omega Model with Omega"};
            obj.omegaModelNames = {"Fixed Omega Model", "Dynamic Omega Model with Qoppa", "Dynamic Omega Model with Omega"};
            folderPath = "/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/models/OmegaModels";
            fileList = dir(fullfile(folderPath, "*.m"));
            obj.numOmegas = numel(fileList);
            obj.models = {
                Model(obj),
                GoBiasModel(obj),
                FixedMotivationalBiasModel(obj),
                DynamicMotivationalBiasModel(obj),
                FixedOmegaModel(obj),
                DynamicOmegaModel1(obj),
                DynamicOmegaModel2(obj),
            };

        end
    end
end