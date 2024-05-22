function [loglik] = model_goBias(parameters,subj)
nd_ep = parameters(1);
ep = 1/(1+exp(-nd_ep));

nd_rho = parameters(2);
rho = exp(nd_rho);

goBias = parameters(3);

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.states;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;

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

        p1 = 1./(1+exp(-(w_g(s)-w_ng(s))));
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

