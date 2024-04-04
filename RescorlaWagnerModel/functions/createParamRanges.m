function paramRanges = createParamRanges(epsilonRange, betaRange, rhoRange)
    arguments
        epsilonRange = [0.01, 0.9];
        betaRange = [1, 10];
        rhoRange = [0.01, 5];
    end

    epsilons = linspace(epsilonRange(1), epsilonRange(2), 10);
    betas = linspace(betaRange(1), betaRange(2), 10);
    rhos = linspace(rhoRange(1), rhoRange(2), 10);

    paramRanges = {epsilons, betas, rhos};

end