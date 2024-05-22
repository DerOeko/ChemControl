function [loglik] = model_dynamicOmega2(parameters,subj)
nd_ep = parameters(1);
ep = 1/(1+exp(-nd_ep));

nd_rho = parameters(2);
rho = exp(nd_rho);

goBias = parameters(3);

nd_oi = parameters(4);
oi = 1/(1+exp(-nd_oi));

nd_alpha = parameters(5);
alpha = 1/(1+exp(-nd_alpha));

nd_beta = parameters(6);
beta = exp(nd_beta);

nd_thres = parameters(7);
thres = tanh(nd_thres);

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.states;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5] * rho;
omega = oi;
Omega = 0;
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

        w_g(s) = (1-omega) * q_g(s) + goBias + omega * sv(s);
        w_ng(s) = (1-omega) * q_ng(s);

        p1 = 1./(1+exp(-(w_g(s)-w_ng(s))));
        p2 = 1-p1;
        
        v_pe = rho * o - sv(s);
        sv(s) = sv(s) + ep * (rho * o - sv(s));
        
        if a==1
            loglik = loglik + log(p1 + eps);
            %q_pe = sqrt((rho * o - q_g(s)).^2 + eps);
            q_pe = rho*o-q_g(s);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
        elseif a==2
            loglik = loglik + log(p2 + eps);
            %q_pe = sqrt((rho * o - q_ng(s)).^2 + eps);
            q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
        end

        Omega = Omega + alpha*(q_pe - v_pe - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));
    end
end
end

