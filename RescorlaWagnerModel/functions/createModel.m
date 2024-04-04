function model = createModel(Qinit, Vinit, goBias, V, pi, modelType, epsilon, beta, rho)
    arguments
        Qinit
        Vinit
        goBias
        V
        pi
        modelType = "Basic"
        epsilon = 0.1
        beta = 1
        rho = 1
    end

    if strcmp(modelType, "Go Bias")
        model = GoBiasModel(epsilon, rho, beta, Qinit, goBias);
    elseif strcmp(modelType, "Fixed Motivational Bias")
        model = FixedMotivationalBiasModel(epsilon, rho, beta, Qinit, goBias, V, pi);
    elseif strcmp(modelType, "Dynamic Motivational Bias")
        model = DynamicMotivationalBiasModel(epsilon, rho, beta, Qinit, goBias, Vinit, pi);
    else
        model = Model(epsilon, rho, beta, Qinit);
    end
end


        
    

