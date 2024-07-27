function [loglik] = ChemControl_mod23(parameters,subj)
% M03 + average reward rate

% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
pi = parameters(4);
alpha_lr = sigmoid(parameters(5));
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
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5] * rho;

loglik = 0;
mu = 0;
for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    sv = initV;

    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        mu = mu + alpha_lr * (o - mu);
        w_g(s) = q_g(s) + goBias + pi * sv(s);
        w_ng(s) = q_ng(s);

        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = rho * o - q_g(s) + mu;
            q_g(s) = q_g(s) + ep * q_pe;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_pe = rho * o - q_ng(s) + mu;
            q_ng(s) = q_ng(s) + ep * q_pe;
        end
    end
end
end

