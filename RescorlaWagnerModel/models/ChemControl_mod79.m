function [loglik] = ChemControl_mod79(parameters,subj)

% Proposal Romain 09/10
% Same as target with perseveration/alternation parameter

% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
epsilon = sigmoid(parameters(4));
betaControl = exp(parameters(5));
betaOmega = exp(parameters(6));
betaPersist = parameters(7);

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

persist_flip=[1,-1];

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

    for t=1:T

        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);

        if t==1
            w_g(s) = q_g(s) + goBias + sv(s)/(betaControl+betaOmega*Omega);
            w_ng(s) = q_ng(s) + (-sv(s))/(betaControl+betaOmega*Omega);
        else
            w_g(s) = q_g(s) + goBias + (sv(s)/(betaControl+betaOmega*Omega))+persist_flip(actions(b, t-1))*betaPersist;
            w_ng(s) = q_ng(s) + ((-sv(s))/(betaControl+betaOmega*Omega))-persist_flip(actions(b, t-1))*betaPersist;
        end
            
        p1 = (epsilon/2) + (1-epsilon)*stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        v_pe = o - sv(s)/rho;

        sv(s) = sv(s) + ep * v_pe;
        
        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = o-q_g(s)/rho;
            q_g(s) = q_g(s) + ep * q_pe;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_pe = o-q_ng(s)/rho;
            q_ng(s) = q_ng(s) + ep * q_pe;
        end

        Omega = Omega + ep*(abs(v_pe)-abs(q_pe) - Omega);
    end
end
end

