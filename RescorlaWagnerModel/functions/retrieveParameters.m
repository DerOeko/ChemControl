function optimalParams = retrieveParameters(modelType, paramRanges)
    arguments
        modelType = "Basic";
        % Parameter ranges
        paramRanges = dictionary("epsilons", 0.01:0.02:0.8, "betas", 1:1:15, "rhos", 0.01:1:15)
    end