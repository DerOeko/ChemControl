 function [loglik] = ChemControl_mod23(parameters,subj)

% Standard Rescorla Wagner model with rho feedback sensitivity + goBias +
% fixed Pavlov + average reward rate modulation + scaled by rho + influence
% on learning rate
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

% Number of blocks:
B = size(outcomes, 1);

% Number of trials:
T = size(outcomes, 2);
initQ = [0 0 0 0];

loglik = 0;
mu = 0;
% Store actions, outcomes and stimuli

% ----------------------------------------------------------------------- %
%% Calculating log likelihood for action sequence with this model:

for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    sv = [0.5 -0.5 0.5 -0.5];

    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        mu = mu + alpha_lr*(rho*o-mu); % -> 1*rho if optimal behaviour
% Safeguard ep calculation
        ep_new = ep * (1 + (mu / (rho + eps))); % Add eps to rho to prevent division by zero
        % Could be negative? Does not make sense
        % Ensure ep stays within a reasonable range (optional, depending on model)
        ep_new = max(min(ep_new, 1-eps), eps); % Keep ep within (eps, 1-eps)
        ep = ep_new;        
        w_g(s) = q_g(s) + goBias + pi * sv(s);
        w_ng(s) = q_ng(s);
        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;
        
        if a==1
            loglik = loglik + log(p1 + eps);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s) + mu);
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s) + mu);
        end
    end
end
end

