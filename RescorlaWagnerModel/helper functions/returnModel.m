function model = returnModel(model_name, parameters)
    %RETURNMODEL Returns a model based on the model_name specified and the
    %parameters. Useful for parameter retrieval to get the model with the
    %best fitting parameters.
    run('config.m');
    switch model_name
        case 'Model'
            % Model(epsilon, rho)
            model = Model(parameters(1), parameters(2));
        
        case 'GoBiasModel'
            % GoBiasModel(epsilon, rho, Qinit, goBias)
            model = GoBiasModel(parameters(1), parameters(2), Qinit, parameters(3));
        
        case 'FixedPavlovModel'
            % FixedMotivationalBiasModel(epsilon, rho, Qinit, goBias, V, pi)
            model = FixedMotivationalBiasModel(parameters(1), parameters(2), Qinit, parameters(3), V, parameters(4));
        
        case 'DynamicPavlovModel'
            % DynamicMotivationalBiasModel(epsilon, rho, Qinit, goBias, Vinit, pi)
            model = DynamicMotivationalBiasModel(parameters(1), parameters(2), Qinit, parameters(3), Vinit, parameters(4));
        
        case 'FixedOmegaModel'
            % FixedOmegaModel(epsilon, rho, Qinit, goBias, Vinit, omega)
            model = FixedOmegaModel(parameters(1), parameters(2), Qinit, parameters(3), Vinit, parameters(4));
        
        case 'DynamicOmega1Model'
            % DynamicOmegaModel1(epsilon, rho, Qinit, goBias, Vinit, omegaInit, alpha, kappa)
            model = DynamicOmegaModel1(parameters(1), parameters(2), Qinit, parameters(3), Vinit, parameters(4), parameters(5), parameters(6));
        
        case 'DynamicOmega2Model'
            % DynamicOmegaModel2(epsilon, rho, Qinit, goBias, Vinit, omegaInit, OmegaInit, alphaOmega, betaOmega, thresOmega)
            model = DynamicOmegaModel2(parameters(1), parameters(2), Qinit, parameters(3), Vinit, parameters(4), OmegaInit, parameters(5), parameters(6), parameters(7));
        
        otherwise
            error('Model not found.');
    end
end

