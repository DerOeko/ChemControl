function [loglik] = ChemControl_mod7(parameters,subj)

% Standard Q learning model with rho feedback sensitivity
% + GoBias
% + DynamicPavlovianBias
% + Initial Optimism
% + ...
% ----------------------------------------------------------------------- %
%% Retrieve parameters:
ep = sigmoid(parameters(1));
rho = exp(parameters(2));
goBias = parameters(3);
omegaBias = -2+4*sigmoid(parameters(4)); % constant bias in control inference
initialOptimism=sigmoid(parameters(5));
% ----------------------------------------------------------------------- %
epsilon = 0;

%% Bayesian module
% initialization of learning (all this can be turned in free parameters)
init_msas = 0.5;   % prior mean, instrumental
init_psas = 1;   % prior confidence, instrumental
init_mss = 0.5;   % prior mean, Pavlovian
init_pss = 1;   % prior confidence, Pavlovian
omegaStartingBias=0; % starting bias in control inference

% switch to activate counterfactual learning
omega_counterfactual=true;

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


favorable = subj.outcomes;
favorable(isNonWinState & favorable == 0) = 1;
favorable(isNonWinState & favorable == -1) = 0;

% Number of blocks:
B = size(outcomes, 1);

% Number of trials:
T = size(outcomes, 2);
initQ = [0 0 0 0];

loglik = 0;

% initializations for MB-Bayes decision making
% for probabilities of favorable outcome
initSAS = zeros(4,2)+init_msas;
initSS = zeros(4,1)+init_mss;
% for confidences into the former
initSASp = init_psas+zeros(4,2);
initSSp = init_pss+zeros(4,1);
% for potential outcomes of gain and loss trials assumed to be known by MB
% Bayes decision-maker (not something to be updated)
initRval = zeros(4,2);
initRval(:,1)=[0 -1 0 -1]; 
initRval(:,2)=[1 0 1 0];

% ----------------------------------------------------------------------- %
%% Calculating log likelihood for action sequence with this model:

for b = 1:B
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    sv = [initialOptimism 1-initialOptimism initialOptimism 1-initialOptimism];

    mss=initSS;
    msas=initSAS;
    
    Mss=initSSp;
    Msas=initSASp;

    Lglob= omegaStartingBias;

    for t=1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        fav = favorable(b, t);

        % update omega values
        omegaGlob =1 / (1 + exp(-omegaBias - Lglob));
        
        % model based values: (Bayesian prob)*(known magnitude) for each action and outcome (favorable or not)
        mb_g=(1-msas(s,1))*initRval(s,1)+msas(s,1)*rho*initRval(s,2);
        mb_ng=(1-msas(s,2))*initRval(s,1)+msas(s,2)*rho*initRval(s,2);
        
        w_g(s) = goBias + (q_g(s)+sv(s))*(1-omegaGlob)+ mb_g*omegaGlob;
        w_ng(s) = (-sv(s)+q_ng(s))*(1-omegaGlob)+ mb_ng*omegaGlob;

        p1 = (epsilon/2) + (1-epsilon)*stableSoftmax(w_g(s), w_ng(s));
        p2 = 1-p1;

        sv(s) = sv(s) + ep * (rho * o - sv(s));
        
        if a==1
            loglik = loglik + log(p1 + eps);
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
        end

        % update controllability likelihood (positive <=> control)
        if fav == 1
            Lglob = Lglob + log(msas(s,a)) - log(mss(s));
        else
            Lglob = Lglob + log(1-msas(s,a)) - log(1-mss(s));
        end

        % update Bayesian transition probabilities
        Mss(s) = Mss(s) + 1;
        Msas(s,a) = Msas(s,a) + 1;
        mss(s) = mss(s) + (fav-mss(s))/Mss(s);
        msas(s,a) = msas(s,a) + (fav-msas(s,a))/Msas(s,a);
        
        if omega_counterfactual
            % Current action a is either 1 or 2, so we calculate the counterfactual action
            counterfactual_action = 3 - a;
            
            % Update counterfactual action msas for state s
            msas(s, counterfactual_action) = omegaGlob * (1 - msas(s, a)) + (1 - omegaGlob) * msas(s, counterfactual_action);
            
            % Update Msas for the counterfactual action
            Msas(s, counterfactual_action) = Msas(s, counterfactual_action) + omegaGlob;
        end
    end
end
end

