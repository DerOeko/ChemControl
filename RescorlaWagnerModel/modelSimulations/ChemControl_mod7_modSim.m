function [out] = ChemControl_mod7_modSim(parameters, subj)
% Standard Q-learning model with delta learning rule.
    
    % ----------------------------------------------------------------------- %
    %% Retrieve parameters:
    ep = sigmoid(parameters(1));
    rho = exp(parameters(2));
    gB = parameters(3);
    alpha = sigmoid(parameters(4));
    beta = exp(parameters(5));
    thres = tanh(0.5*parameters(6));

    % ----------------------------------------------------------------------- %
    %% Unpack data:

    % Extract task features:
    stimuli = subj.stimuli; 
    controllabilities = subj.controllability; % 1 or 0
    randomRewards = subj.randomReward; % 1 or 0
    randHCs = subj.randHC; % 1, 0, 2
    randLCs = subj.randLC; % 1, 0, 2
    cali_stimuli = subj.cali_stimuli;
    cali_randHC = subj.cali_randHC;
    cali_randLC = subj.cali_randLC;
    cali_randRewards = subj.cali_randReward;


    % Data dimensions:
    B = size(stimuli, 1); % Number of blocks
    T = size(stimuli, 2); % Number of trials
    S = 4; % Number of stimuli

    % Store outputs
    HCcell = cell(B/2, S); % Store go probs in HC
    LCcell = cell(B/2, S); % Store go probs in LC
    actions = zeros(B, T);
    outcomes = zeros(B, T);

    q0 = [0.5 -0.5 0.5 -0.5];
    hc = 0;
    lc = 0;
    

    %% Run calibration block
    rewardLossCounter = zeros([1, 2]);
    q_g = q0 * rho;
    q_ng = q0 * rho;
    w_g = q0 * rho;
    w_ng = q0 * rho;
    sv = q0;

    isHC = 1;
    Omega = 0;
    omega = 1/(1+exp(-beta*(Omega-thres)));

    for t = 1:T
        s = cali_stimuli(t);
        randHC = cali_randHC(t);
        randLC = cali_randLC(t);
        isRewarded = cali_randRewards(t);

        w_g(s) = (1-omega) * q_g(s) + gB + omega * sv(s);
        w_ng(s) = (1-omega) * q_ng(s);

        p1 = stableSoftmax(w_g(s), w_ng(s));

        a = returnAction(p1);
        o = returnReward(s, a, isHC, randLC, randHC, isRewarded);

        v_pe = rho * o - sv(s);
        sv(s) = sv(s) + ep * (rho * o - sv(s));

        if a==1
            q_pe = rho*o-q_g(s);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
        elseif a==2
            q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
        end

        Omega = Omega + alpha*(q_pe - v_pe - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));

        counter = updateRewardLossCounter(s, o);
        rewardLossCounter = rewardLossCounter + counter;
    end
   
    M = rewardLossCounter/(T/2);
    averageRewardRate = M(1);
    averageLossRate = M(2);
    numRewarded = round(averageRewardRate * T); % number of rewarded trials
    numAvoided = round(averageLossRate * T);

    %% Simulating data
    for b = 1:B
        switch controllabilities(b, 1)
            case 1
                isHC = true;
                isLC = false;
                isYoked = false;
            case 0
                isHC = false;
                isLC = true;
                isYoked = false;
            case 2
                isHC = false;
                isLC = true;
                isYoked = true;

                rewardedVec = [ones(1, numRewarded) zeros(1, T-numRewarded, 1)];
                rewardedVec = rewardedVec(randperm(length(rewardedVec)));
                avoidedVec = [ones(1, numAvoided) zeros(1, T-numAvoided, 1)];
                avoidedVec = avoidedVec(randperm(length(avoidedVec)));
        end
        hc = hc + isHC;
        lc = lc + isLC;
        
        q_g = q0 * rho;
        q_ng = q0 * rho;
        w_g = q0 * rho;
        w_ng = q0 * rho;
        sv = q0;

        Omega = 0;
        omega = 1/(1+exp(-beta*(Omega-thres)));

        for t = 1:T
            s = stimuli(b, t);
            isWinState = mod(s, 2);
            randHC = randHCs(b, t); % outcome matters (1, 0, 2)
            randLC = randLCs(b, t); % outcome doesn't matter (1, 0, 2)
            if ~isYoked
                isRewarded = randomRewards(b, t);
            elseif isWinState && isYoked
                isRewarded = rewardedVec(t);
            elseif ~isWinState && isYoked
                isRewarded = avoidedVec(t);
            end

            w_g(s) = (1-omega) * q_g(s) + gB + omega * sv(s);
            w_ng(s) = (1-omega) * q_ng(s);


        p1 = stableSoftmax(w_g(s), w_ng(s));
            
            if isHC
                HCcell{hc, s}(end+1) = p1;
            else
                LCcell{lc, s}(end+1) = p1;
            end

            a = returnAction(p1);
            o = returnReward(s, a, isHC, randLC, randHC, isRewarded);
            actions(b, t) = a;
            outcomes(b, t) = o;

            v_pe = rho * o - sv(s);
            sv(s) = sv(s) + ep * (rho * o - sv(s));
            if a==1
                q_pe = rho*o-q_g(s);

                q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
            elseif a==2
                q_pe = rho*o-q_ng(s);

                q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
            end

            Omega = Omega + alpha*(q_pe - v_pe - Omega);
            omega = 1/(1+exp(-beta*(Omega-thres)));
        end
    end
    
    out.HCcell = HCcell;
    out.LCcell = LCcell;
    out.randHC = randHCs;
    out.randLC = randLCs;
    out.stimuli = stimuli;
    out.controllability = controllabilities;
    out.randomRewards = randomRewards;
    out.actions = actions;
    out.outcomes = outcomes;
end