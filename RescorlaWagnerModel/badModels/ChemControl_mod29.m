function [loglik] = ChemControl_mod29(parameters,subj)

% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
% ----------------------------------------------------------------------- %

%% Unpack data:
actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

% Number of blocks:
B = size(outcomes, 1);

% Number of trials:
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5];

loglik = 0;

% Store actions, outcomes and stimuli

% ----------------------------------------------------------------------- %
%% Calculating log likelihood for action sequence with this model:

for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;

    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);

        w_g(s) = q_g(s);
        w_ng(s) = q_ng(s);
        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        if a==1
            loglik = loglik + log(p1 + eps);
            q_g(s) = q_g(s) + ep * (o - q_g(s));
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_ng(s) = q_ng(s) + ep * (o - q_ng(s));
        end
    end
end
end

