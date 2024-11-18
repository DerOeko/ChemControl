function [loglik] = ChemControl_mod2(parameters,subj)
% Standard Q learning model with rho feedback sensitivity
% + GoBias
% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
% ----------------------------------------------------------------------- %
%% Unpack data:

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;
% Identify win and non-win states
isWinState = mod(states, 2) == 1;
isNonWinState = ~isWinState;

% Transform outcomes for win states
outcomes(isWinState & outcomes == 0) = -1;

% Transform outcomes for non-win states
outcomes(isNonWinState & outcomes == 0) = 1;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0 0 0 0];

loglik = 0;

for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;

    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);

        w_g(s) = q_g(s) + goBias;
        w_ng(s) = q_ng(s);

        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        if a==1
            loglik = loglik + log(p1 + eps);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
        end
    end
end
end

