function [loglik] = ChemControl_mod4(parameters, subj)
% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
pi = parameters(4);
% ----------------------------------------------------------------------- %
%% Unpack data:
actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5];

loglik = 0;

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

        w_g(s) = q_g(s) + goBias + pi * sv(s);
        w_ng(s) = q_ng(s);

        p1 = stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;
        
        sv(s) = sv(s) + ep * (rho * o - sv(s));
        
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

