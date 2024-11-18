function [loglik] = ChemControl_mod87(parameters, subj)
% Model, inspired by Sam Gershman paper

% Parameter Transformation
bt = exp(parameters(1));         % Inverse temperature
mq = sigmoid(parameters(2));     % Prior mean instrumental
pq = exp(parameters(3));        % Prior confidence instrumental
mv = sigmoid(parameters(4));     % Prior mean Pavlovian
pv = exp(parameters(5));        % Prior confidence Pavlovian
w0 = sigmoid(parameters(6));     % Initial Pavlovian weight (log odds)

% Unpack data
actions = subj.actions;
outcomes = subj.outcomes;
states = subj.stimuli;

% Preprocess outcomes
isWinState = mod(states, 2) == 1;       % Identify Win States
isLossState = ~isWinState;             % Identify Loss States

% For Loss States:
% - Change outcome == -1 to 0
% - Change outcome == 0 to 1
outcomes(isLossState & outcomes == -1) = 0;
outcomes(isLossState & outcomes == 0)  = 1;

% For Win States:
% - Outcomes remain as is (0 or 1)
% No changes needed

B = size(outcomes, 1); % Number of blocks
T = size(outcomes, 2); % Number of trials
unique_states = unique(states);
S = length(unique_states);
Q = zeros(S, 2) + mq;          % Q(s,1) and Q(s,2)
V = zeros(S, 1) + mv;           % V(s)
Mq = ones(S, 2) * pq;
Mv = ones(S, 1) * pv;
loglik = 0;

for b = 1:B
    % Initialize log-odds at the start of each block
    L = log(w0 + eps) - log(1 - w0 + eps); % can be reset every block or only once at the beginning

    
    for t = 1:T
        a = actions(b, t);
        o = outcomes(b, t);
        s = states(b, t);
        stim = find(unique_states == s);
        
        
        w = 1 ./ (1 + exp(-L));
        
        d = (1 - w) * Q(stim, 1) - (1 - w) * Q(stim, 2) - w * V(stim);
        
        P = 1 ./ (1 + exp(-bt * d));
        
        if a == 1  % Go
            loglik = loglik + log(P + eps);
        else      % NoGo
            loglik = loglik + log(1 - P + eps);
        end
        
        if o == 1
           
            L = L + log(V(stim) + eps) - log(Q(stim, a) + eps);
        else
          
            L = L + log(1 - V(stim) + eps) - log(1 - Q(stim, a) + eps);
        end
        
        Mv(stim) = Mv(stim) + 1;
        Mq(stim, a) = Mq(stim, a) + 1;

        V(stim) = V(stim) + (o - V(stim)) / Mv(stim);
        Q(stim, a) = Q(stim, a) + (o - Q(stim, a)) / Mq(stim, a);
        
    end
end
end

