function [loglik] = ChemControl_mod20(parameters,subj)

ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
alpha = sigmoid(parameters(4));
beta = exp(parameters(5));
thres = scaledSigmoid(parameters(6));

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5] * rho;

loglik = 0;

for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    sv = initV;
    Omega = 0;
    omega = 1/(1+exp(-beta*(-thres)));
    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);

        w_g(s) = q_g(s) + goBias + sv(s);
        w_ng(s) = q_ng(s);

        p1 = stableSoftmax(w_g(s), w_ng(s));
   
        p2 = 1-p1;
        
        v_pe = o - sv(s);
        ep = ep * (1-omega);
        
        if a==1
            loglik = loglik + log(p1 + eps);
            %q_pe = sqrt((rho * o - q_g(s)).^2 + eps);
            q_pe = o-q_g(s);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
        elseif a==2
            loglik = loglik + log(p2 + eps);
            %q_pe = sqrt((rho * o - q_ng(s)).^2 + eps);
            q_pe = o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
        end

        Omega = Omega + alpha*(abs(q_pe)- abs(v_pe) - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));
    end
end
end

