function [loglik] = ChemControl_mod32(parameters,subj)

% Dynamic Omega, dynamic pavlov, non competitive, ARR on V and
% Q, rescaled PEs
% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
alpha = sigmoid(parameters(4));
beta = exp(parameters(5));
thres = scaledSigmoid(parameters(6));
alpha_lr = sigmoid(parameters(7));

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

        w_g(s) = q_g(s) + goBias + (1-omega) * sv(s);
        w_ng(s) = q_ng(s);
        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        v_pe = o - sv(s)/rho + mu;
        sv(s) = sv(s) + ep * (rho * o - sv(s) + mu);

        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = o-q_g(s)/rho + mu;

            q_g(s) = q_g(s) + ep * (rho * o - q_g(s) + mu);
        elseif a==2
            loglik = loglik + log(p2 + eps);

            q_pe = o-q_ng(s)/rho + mu;
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s) + mu);
        end

        Omega = Omega + alpha*(v_pe - q_pe - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));
    end
end
end

