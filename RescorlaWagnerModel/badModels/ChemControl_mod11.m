function [loglik] = ChemControl_mod11(parameters,subj)

% personal guess model (what I would  plemented to test our
% hypothesis)

ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
alpha_up = ep;
alpha_down = ep; 

actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

B = size(outcomes, 1);
T = size(outcomes, 2);
initQ = [0.5 -0.5 0.5 -0.5] * rho;
initV = [0.5 -0.5 0.5 -0.5] * rho;
initRR=0;

loglik = 0;

for b = 1:B
    
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    rr=initRR;
    sv = initV;
    omega = 0.5;

    for t=1:T
        
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);

        % compute decision weights for go and no go (include the pavlovian
        % values in the no go side of things
        w_g(s) = omega * q_g(s) + (1-omega) * sv(s) + goBias;
        w_ng(s) = omega * q_ng(s) + (1-omega) * (-sv(s));
        
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
            q_pe_abs=abs(o - q_ng(s));
        end

        % update Pavlovian model
        %v_pe = rho * o - sv(s);
        sv(s) = sv(s) + ep * (rho * o - sv(s));
        
        % update Instrumental model
        if a==1
            loglik = loglik + log(p1 + eps);
            %q_pe = rho*o-q_g(s);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
            %p_explore=p2;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            %q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
            %p_explore=p1;
        end
        
        % update global reward rate
        %rr=rr+ep*(rho * o - sv(s));
        sigdiff=sigmoid((v_pe_abs-q_pe_abs)*1000);
        % different logics to update Omega: simply tracks the overall frequency of trials
        % where the prediction error of the Pavlovian model is higher than that of the instrumental model
        omega = omega + (sigdiff)*(alpha_up*(1-omega)) + (1-sigdiff) *(alpha_down*(0-omega));
        % if v_pe_abs>q_pe_abs
        %     omega = omega + alpha_up*(1 - omega);
        % elseif v_pe_abs<q_pe_abs
        %     omega = omega + alpha_down*(0 - omega);
        % else
        %     omega = omega;
        % end     
    end
end
end

