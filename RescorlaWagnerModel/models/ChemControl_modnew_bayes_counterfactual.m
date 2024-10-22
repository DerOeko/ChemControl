function [loglik] = ChemControl_mod13(parameters,subj)

% Proposal Romain 09/10
% This is just the same models as the target model, but with the dynamic omega effect removed

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
omega_counterfactual=true

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

%%% make favorable outcome vector TO DOUBLE CHECK
% should be 1 for any favorable outcome and 0 for any unfavorable outcome
% not sure about the logic of your input outcome vector so this is a guess
favorable=subj.outcomes;
favorable(isNonWinState & favorable == 0) = 1;
favorable(isWinState & favorable == 0) = 0;

% Number of blocks:
B = size(outcomes, 1);

% Number of trials:
T = size(outcomes, 2);
initQ = zeros(4,1);
initV=[initialOptimism 1-initialOptimism initialOptimism 1-initialOptimism];

% initializations for MB-Bayes decision making
%
initSAS = zeros(4,2)+init_msas;
initSS = zeros(4,1)+init_mss;
%
initSASp = init_psas+zeros(4,2);
initSSp = init_pss+zeros(4,1);
% 
initRval = zeros(4,2);
initRval(:,1)=[0 -1 0 -1]; % assumes that participants know potential outcomes of gain and loss trials
initRval(:,2)=[1 0 1 0];
%

loglik = 0;

% Store actions, outcomes and stimuli

% ----------------------------------------------------------------------- %
%% Calculating log likelihood for action sequence with this model:
Omega = 0;

for b = 1:B
    
    w_g = initQ;
    w_ng = initQ;
    q_g = initQ;
    q_ng = initQ;
    sv = initV;
    
    mss=initSS;
    msas=initSAS;
    
    Mss=initSSp;
    Msas=initSASp;
    
    Lglob=omegaStartingBias;
    Lstate=[0 0 0 0]+omegaStartingBias;
    omegaState=[0 0 0 0];

    for t=1:T

        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        fav = favorable(b, t);
        
        % update omega values
        omegaGlob =1 / (1 + exp(-omegaBias - Lglob));
        for st=1:4
            omegaState(s) =1 / (1 + exp(-omegaBias - Lstate(s))); % we could instead just compute current omegaState(s), not all 4
        end

        % model based values: (Bayesian prob)*(known magnitude) for each action and outcome (favorable or not)
        mb_g=msas(s,1,1)*initRval(s,1)+msas(s,1,2)*rho*initRval(s,2);
        mb_ng=msas(s,2,1)*initRval(s,1)+msas(s,2,2)*rho*initRval(s,2);
        
        w_g(s) = goBias + (q_g(s)+sv(s))*(1-omega(s))+ mb_g*omega(s);
        w_ng(s) = (-sv(s)+q_ng(s))*(1-omega(s))+ mb_ng*omega(s);
        
        p1 = (epsilon/2) + (1-epsilon)*stableSoftmax(w_g(s), w_ng(s));

        p2 = 1-p1;

        v_pe = rho*o - sv(s);
        sv(s) = sv(s) + ep * v_pe;

        if a==1
            loglik = loglik + log(p1 + eps);
            q_pe = rho*o-q_g(s);
            q_g(s) = q_g(s) + ep * q_pe;
        elseif a==2
            loglik = loglik + log(p2 + eps);
            q_pe = rho*o-q_ng(s);
            q_ng(s) = q_ng(s) + ep * q_pe;
        end
        
        
        % update model posterior
        if fav == 1
            Lglob = Lglob + log(msas(s,a)) - log(mss(s));
            Lstate(s)=Lstate(s) + log(msas(s,a)) - log(mss(s));
        else
            Lglob = Lglob + log(1-msas(s,a)) - log(1-mss(s));
            Lstate(s)=Lstate(s) + log(1-msas(s,a)) - log(1-mss(s));
        end
                
        % update Bayesian transition probabilities
        Mss(s) = Mss(s) + 1;
        Msas(s,a) = Msas(s,a) + 1;
        mss(s) = mss(s) + (fav-mss(s))/Mv(s);
        msas(s,a) = msas(s,a) + (fav-msas(s,a))/Mq(s,a);
        
        if omega_counterfactual
            msas(s,2-a) = omegaGlob*(1-msas(s,a)) + (1-omegaGlob)*msas(s,2-a);
            Msas(s,2-a) = Msas(s,2-a) + omegaGlob;
        end

    end
end
end

