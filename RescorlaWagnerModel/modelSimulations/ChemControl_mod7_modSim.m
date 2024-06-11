function [out] = ChemControl_mod7_modSim(parameters, subj)
% Standard Q-learning model with delta learning rule.
    
    % ----------------------------------------------------------------------- %
    %% Retrieve parameters:
    ep = sigmoid(parameters(1));
    rho = exp(parameters(2));
    gB = parameters(3);
    oi = sigmoid(parameters(4));
    alpha = sigmoid(parameters(5));
    beta = exp(parameters(6));
    thres = tanh(parameters(7));

    % ----------------------------------------------------------------------- %
    %% Unpack data:

    % Extract task features:
    stimuli = subj.stimuli; 
    controllabilities = subj.controllability; % 1 or 0
    randomRewards = subj.randomReward; % 1 or 0
    randHCs = subj.randHC; % 1, 0, 2
    randLCs = subj.randLC; % 1, 0, 2

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
    
    omega = oi;
    Omega = 0;
    %% Simulating data
    for b = 1:B
        isHC = controllabilities(b, 1);
        hc = hc + isHC;
        lc = lc + ~isHC;
        
        q_g = q0 * rho;
        q_ng = q0 * rho;
        w_g = q0 * rho;
        w_ng = q0 * rho;
        sv = q0;

        for t = 1:T
            s = stimuli(b, t);
            randHC = randHCs(b, t); % outcome matters (1, 0, 2)
            randLC = randLCs(b, t); % outcome doesn't matter (1, 0, 2)
            isRewarded = randomRewards(b, t);

            w_g(s) = (1-omega) * q_g(s) + gB + omega * sv(s);
            w_ng(s) = (1-omega) * q_ng(s);


            p1 = exp(w_g(s))./(exp(w_g(s)) + exp(w_ng(s)));
            
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