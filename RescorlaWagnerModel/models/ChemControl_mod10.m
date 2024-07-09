function [loglik] = ChemControl_mod10(parameters,subj)

% personal guess model (what I would have implemented to test our
% hypothesis)

ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
alpha = sigmoid(parameters(4));
beta = exp(parameters(5));
thres = scaledSigmoid(parameters(6));

% additional parameters
w_rew_info=parameters(7); % may go both ways
% omega_bounds=0.5*sigmoid(parameters(8)); % cannot be bigger than 0.5
% softmax_bounds=0.5*sigmoid(parameters(8)); % cannot be bigger than 0.5


actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5];
initRR=0;

loglik = 0;

for b = 1:B
    
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    rr=initRR;
    sv = initV;
    Omega = 0;

    for t=1:T
        
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        
        % omega is assessed before each decision, in a state dependent
        % fashion
        % note that higher Omega reflect higher control
        % the w_rew_info reflects the relative importance of expected reward versus
        % information in the estimation of controllability
        % reward expectations are renormalized before -1 to 1
        Xomega=-beta*(w_rew_info*(sv(s)/rho)+(1-w_rew_info)*Omega-thres);
        
        % to test a version that biases omega as a function of reward rate,
        % use instead:
        % Xomega=-beta*(w_rew_info*(rr/rho)+(1-w_rew_info)*Omega-thres);
        
        % create the small omega arbitrator
        omega=1/(1+exp(Xomega));
        
        % to test a version that bounds the arbitrator between [omega_bounds, 1-omega_bounds] 
        % use instead:
        % omega = omega_bounds+((1-2*omega_bounds)/(1+exp(Xomega)));

        % compute decision weights for go and no go (include the pavlovian
        % values in the no go side of things
        w_g(s) = omega * q_g(s) + (1-omega) * sv(s) + goBias;
        w_ng(s) = omega * q_ng(s) + (1-omega) * (-sv(s));
        
        % alternative weighting 1: the go bias is expressed more when omega
        % is low
        % w_g(s) = omega * q_g(s) + (1-omega) * (sv(s) + goBias);
        % w_ng(s) = omega * q_ng(s) + (1-omega) * (-sv(s));
        
        % alternative weighting 2: the go bias becomes a general pavlovian
        % bias scaling how much Pavlovian values impact go vs no_go
        % w_g(s) = omega * q_g(s) + (1-omega) * (sv(s)*goBias);
        % w_ng(s) = omega * q_ng(s) + (1-omega) * (-sv(s)*goBias);

        % compute decision probabilities
        p1 = stableSoftmax(w_g(s), w_ng(s));
    
        % to implement random exploration, use instead:
        % p1 = stableSoftmax_randexp(w_g(s), w_ng(s),softmax_bounds);
        
        % no go prob
        p2 = 1-p1;
        
        % compute absolute reward prediction errors (between 0 and 1)
        v_pe_abs=abs(o - sv(s));
        if a==1
            q_pe_abs=abs(o - q_g(s));
        else
            q_pe_abs=abs(o-q_ng(s));
        
        % update Pavlovian model
        v_pe = rho * o - sv(s);
        sv(s) = sv(s) + ep * (rho * o - sv(s));
        
        % update Instrumental model
        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = rho*o-q_g(s);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
            p_explore=p2;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
            p_explore=p1;
        end
        
        % update global reward rate
        rr=rr+ep*(rho * o - sv(s));

        % update Omega in a way that scales with the "probability" that the
        % participant has made an exploratory decision on that trial
        Omega = Omega + (alpha*p_explore)*(abs(v_pe)-abs(q_pe) - Omega);
        
        % to update continuously, use instead
        % Omega = Omega + alpha*(abs(v_pe)-abs(q_pe) - Omega);
                
    end
end
end

