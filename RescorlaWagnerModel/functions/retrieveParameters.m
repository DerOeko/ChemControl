function optimalParams = retrieveParameters(simulatedActions, stimulusSequence, outcomeVecs, successVecs, correctActions, modelType, paramRanges)
    arguments
        simulatedActions
        stimulusSequence
        outcomeVecs
        successVecs
        correctActions
        modelType = "Basic";
        % Parameter ranges
        paramRanges = createParamRanges()
    end
    
    epsilonRange = paramRanges{1};
    betaRange = paramRanges{2};
    rhoRange = paramRanges{3};

    bestParams = [];
    bestFit = Inf;
    
    for epsilon = epsilonRange
        for beta = betaRange
            for rho = rhoRange
                model = createModel(Qinit, Vinit, goBias, V, pi, modelType, epsilon, beta, rho);
                negLogLikelihood = calculateFit(model, simulatedActions, stimulusSequence, outcomeVecs, successVecs, correctActions);

                if negLogLikelihood < bestFit
                    bestParams = [epsilon, beta, rho];
                    bestFit = negLogLikelihood;
                end
            end
        end
    end
    optimalParams = bestParams;

