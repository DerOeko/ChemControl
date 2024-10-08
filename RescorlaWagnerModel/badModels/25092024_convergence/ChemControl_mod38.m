function [loglik] = ChemControl_mod38(parameters,subj)

% Dynamic Omega, dynamic pavlov, without pavlov on the nogo side,
% p_explore, arr, randexp by omega, with ep = alpha_lr = alpha, fixed rho
% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = 1;
goBias = parameters(2);
alpha = ep;
beta = exp(parameters(3));
thres = scaledSigmoid(parameters(4));
alpha_lr = ep;
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

% Number of blocks:
B = size(outcomes, 1);

% Number of trials:
T = size(outcomes, 2);
initQ = [0 0 0 0];

loglik = 0;

% Store actions, outcomes and stimuli

% ----------------------------------------------------------------------- %
%% Calculating log likelihood for action sequence with this model:
Omega = 0;
mu = 0;
for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    sv = [0.5 -0.5 0.5 -0.5];
    omega = 1/(1+exp(-beta*(Omega-thres)));
    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        mu = mu + alpha_lr*(rho*o-mu);
        w_g(s) = omega * q_g(s) + goBias + (1-omega) * sv(s);
        w_ng(s) = omega * q_ng(s);
        p1 = stableSoftmax_randexp(w_g(s), w_ng(s), 1-omega);
        p2 = 1-p1;

        v_pe = rho*o - sv(s)+mu;
        sv(s) = sv(s) + ep * (rho * o - sv(s)+mu);

        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = rho*o-q_g(s)+mu;

            q_g(s) = q_g(s) + ep * (rho * o - q_g(s)+mu);
            p_explore = p2;
        elseif a==2
            loglik = loglik + log(p2 + eps);

            q_pe = rho*o-q_ng(s)+mu;
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s)+mu);
            p_explore = p1;
        end

        Omega = Omega + (alpha*p_explore)*(v_pe - q_pe - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));
    end
end
end

