function [loglik] = ChemControl_mod14(parameters,subj)
% Model idea by Samu
% Controllability arbitration with average reward modulated prediction
% error signal.

ep = sigmoid(parameters(1));
softmax_beta = exp(parameters(2));
goBias = parameters(3);
alpha = sigmoid(parameters(4));
beta = exp(parameters(5));
thres = scaledSigmoid(parameters(6));
alpha_lr = sigmoid(parameters(7));

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5];
initV = [0.5 -0.5 0.5 -0.5];

loglik = 0;

mu = 0;
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
        mu = mu + alpha_lr * (o - mu);

        w_g(s) = (1-omega) * q_g(s) + goBias + omega * sv(s);
        w_ng(s) = (1-omega) * q_ng(s);

        p1 = stableSoftmax_with_b(w_g(s), w_ng(s), softmax_beta);
   
        p2 = 1-p1;
        
        v_pe = o - sv(s) + mu;

        sv(s) = sv(s) + ep * v_pe;
        
        if a==1
            loglik = loglik + log(p1 + eps);
            %q_pe = sqrt((rho * o - q_g(s)).^2 + eps);
            q_pe = o - q_g(s) + mu;
            q_g(s) = q_g(s) + ep * q_pe;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            %q_pe = sqrt((rho * o - q_ng(s)).^2 + eps);
            q_pe = o - q_ng(s) + mu;
            q_ng(s) = q_ng(s) + ep * q_pe;
        end
        
        Omega = Omega + alpha*(q_pe - v_pe - Omega);
        omega = 1/(1+exp(-beta*(Omega-thres)));
    end
end
end

