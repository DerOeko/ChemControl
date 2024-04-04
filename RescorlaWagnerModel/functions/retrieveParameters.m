function optimalParams = retrieveParameters(simulatedActions, stimulusSequence, outcomeVecs, successVecs, correctActions, controllabilityArray, Qinit, Vinit, goBias, V, pi, modelType, paramRanges)
    arguments
        simulatedActions
        stimulusSequence
        outcomeVecs
        successVecs
        correctActions
        controllabilityArray
        Qinit = zeros(4,2)
        Vinit = [0.5 -0.5 0.5 -0.5] % Default value for V
        goBias = 0.5
        V = [0.5 -0.5 0.5 -0.5] % Default value for V
        pi = 0.5 % Default value for pi
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
                negLogLikelihood = calculateFit(model, simulatedActions, stimulusSequence, outcomeVecs, successVecs, correctActions, controllabilityArray);
                if epsilon == 0.2 && beta == 3 && rho == 0.1
                    disp("-----------")
                    fprintf('Parameter values: epsilon = %.2f, beta = %.2f, rho = %.2f\n', epsilon, beta, rho);
                    fprintf('Negative log-likelihood: %.2f\n', negLogLikelihood);
                    disp("-----------")
                end
                if negLogLikelihood < bestFit
                    bestParams = [epsilon, beta, rho];
                    bestFit = negLogLikelihood;

                    fprintf('Parameter values: epsilon = %.2f, beta = %.2f, rho = %.2f\n', epsilon, beta, rho);
                    fprintf('Negative log-likelihood: %.2f\n', negLogLikelihood);
                end
            end
        end

        
    end
    fprintf('Optimal parameters: epsilon = %.2f, beta = %.2f, rho = %.2f\n', bestParams);

    optimalParams = bestParams;

