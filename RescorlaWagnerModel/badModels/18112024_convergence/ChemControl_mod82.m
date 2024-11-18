function [loglik] = ChemControl_mod82(parameters, subj)
    % Unpack parameters
    ep = sigmoid(parameters(1));
    rho = exp(parameters(2));
    goBias = parameters(3);
    alpha = sigmoid(parameters(4));
    beta = exp(parameters(5));
    thres = scaledSigmoid(parameters(6));

    % Unpack data
    actions = subj.actions;
    outcomes = subj.outcomes;
    states = subj.stimuli;

    % Preprocess outcomes
    isWinState = mod(states, 2) == 1;
    isNonWinState = ~isWinState;
    outcomes(isWinState & outcomes == 0) = -1;
    outcomes(isNonWinState & outcomes == 0) = 1;

    B = size(outcomes, 1); % Number of blocks
    T = size(outcomes, 2); % Number of trials
    initQ = [0 0 0 0];
    loglik = 0;

    Omega = 0;
    omega = 1./(1+exp(-beta*(Omega-thres)));

    Q_actor_g = initQ;
    Q_actor_ng = initQ;

    V_spectator = [0.5 -0.5 0.5 -0.5];

    Q_combined_g = initQ;
    Q_combined_ng = initQ;

    for b = 1:B
        for t = 1:T
            a = actions(b, t);
            o = outcomes(b, t);
            s = states(b, t);
            
            % Generate action weights
            Q_combined_g(s) = omega * Q_actor_g(s) + (1-omega) * V_spectator(s) + goBias;
            Q_combined_ng(s) = omega * Q_actor_ng(s);

            % Determine action
            p_action = stableSoftmax(Q_combined_g(s), Q_combined_ng(s));
            p_noaction = 1-p_action;

            % Update Actor Q-values
            if a==1
                PE_actor = o - Q_actor_g(s)/rho;

                loglik = loglik + log(p_action + eps);
                Q_actor_g(s) = Q_actor_g(s) + ep * (rho * o - Q_actor_g(s));
            else
                PE_actor = o - Q_actor_ng(s)/rho;

                loglik = loglik + log(p_noaction + eps);
                Q_actor_ng(s) = Q_actor_ng(s) + ep * (rho * o - Q_actor_ng(s));
            end

            % Update Spectator Q-values (Pavlovian)
            if isWinState(s)
                P_go = 1;
            else
                P_go = -1;
            end

            PE_spectator = P_go - V_spectator(s);
            V_spectator(s) = V_spectator(s) + ep * (P_go - V_spectator(s));

            Omega = Omega + alpha * (abs(PE_spectator) - abs(PE_actor) - Omega);
            omega = 1/(1+exp(-beta*(Omega-thres)));
        end
    end
end

