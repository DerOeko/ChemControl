function [loglik] = ChemControl_mod11(parameters,subj)

% Dynamic Omega, dynamic pavlov, with scaled PEs, p_explore, competitive
% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
alpha = sigmoid(parameters(4));
beta = exp(parameters(5));
thres = scaledSigmoid(parameters(6));
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

        w_g(s) = omega * q_g(s) + goBias + (1-omega) * sv(s);
        w_ng(s) = omega * q_ng(s);
        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        v_pe = rho*o - sv(s);
        sv(s) = sv(s) + ep * (rho * o - sv(s));

        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = rho*o-q_g(s);

            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
            p_explore = p2;
        elseif a==2
            loglik = loglik + log(p2 + eps);

            q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
            p_explore = p1;
        end

        Omega = Omega + (alpha*p_explore)*(v_pe - q_pe - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));
    end
end
end

