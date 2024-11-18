function [loglik] = ChemControl_mod13(parameters,subj)

% Proposal Romain 09/10
% This is just the same models as the target model, but with the dynamic omega effect removed

% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
epsilon = sigmoid(parameters(4));
betaControl = exp(parameters(5));

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

    omega = omega_init;

    for t=1:T

        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);

        w_g(s) = q_g(s) + goBias + sv(s)/(betaControl);
        w_ng(s) = q_ng(s) +  (-sv(s))/(betaControl);
        
        p1 = (epsilon/2) + (1-epsilon)*stableSoftmax(w_g(s), w_ng(s));

        p2 = 1-p1;

        v_pe = rho*o - sv(s);
        sv(s) = sv(s) + ep * v_pe;

        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = rho*o-q_g(s);
            q_g(s) = q_g(s) + ep * q_pe;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * q_pe;
        end

        Omega = Omega + ep*(abs(v_pe)-abs(q_pe) - Omega);

        omega = Omega

    end
end
end

