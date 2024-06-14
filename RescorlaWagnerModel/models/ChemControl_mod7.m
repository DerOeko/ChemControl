function [loglik] = ChemControl_mod7(parameters,subj)

ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
oi = sigmoid(parameters(4));
alpha = sigmoid(parameters(5));
beta = exp(parameters(6));
thres = scaledSigmoid(parameters(7));

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5];
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

        p1 = 1/(1+exp(w_ng(s)-w_g(s)));
   
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

