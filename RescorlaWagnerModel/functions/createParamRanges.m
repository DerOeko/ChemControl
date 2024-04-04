function paramRanges = createParamRanges(epsilonRange, betaRange, rhoRange)
    arguments
        epsilonRange = [0.01, 0.8];
        betaRange = [1, 15];
        rhoRange = [0.01, 15];
    end

    epsilons = linspace(epsilonRange(1), epsilonRange(2), 20);
    betas = linspace(betaRange(1), betaRange(2), 15);
    rhos = linspace(rhoRange(1), rhoRange(2), 15);

    paramRanges = {epsilons, betas, rhos};

end