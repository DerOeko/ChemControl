function [HCcell, LCcell, cArr, os, cPs] = runExperiment(model, cfg)
    % Function that actually runs the experiment.
    % Receives a model and a config class.
    % The model can be any model.
    % The config class is defined in Config and contains meta parameters.
    T = cfg.T; % Number of trials in a block
    B = cfg.B; % Number of blocks
    rP = cfg.rP; % Reward probability
    cP = cfg.cP; % Controll probability
    cArr = cfg.cArr; % Controllability array

    % Cues and how often they repeat
    cond = [1, 2, 3, 4]; 
    rep = T / 4;

    % Controllability array
    HC_B = B / 2; % number of hc blocks
    HCcell = cell(HC_B, 4); % cell array for storing proportion of go responses
    LCcell = cell(B - HC_B, 4);
    % Whether or not the sequence of controllability blocks is alternating
    % or not (useful for plotting the overage omega across runs.
    if ~contains(class(model), "Omega")
        cArr = zeros(B, 1);
        cArr(1:HC_B) = 1;
        cArr = cArr(randperm(length(cArr)));
    end

    % Logging cells and arrays
    sSeq = cell(B, 1); % Stimulus sequences
    cVecs = cell(B, 1); % Controllability arrays
    rVecs = cell(B, 1); % Reward arrays (whether the current trial is rewarded or not)
    aVec = zeros(B, T); % Action arrays
    correct_aVec = zeros(B, T); % Correct Actions arrays
    os = zeros(B*T, 1); % Omegas arrays
    cPs = zeros(B*T, 1); % Controllability probability arrays
    rMat = zeros(B, T); % Matrix of actual received reward
    
    isO = contains(class(model), "Omega"); % Flag for whether the model is omega
    tIdx = 0; % trial index
    hc = 0; % high control block number
    lc = 0; % low control block number
    m = model; % model

    for b = 1:B % loop over blocks
        % Whether the block is high control
        isHC = cArr(b); % Flag for whether current block is high control 
        
        if isHC
            hc = hc + 1;
        else
            lc = lc + 1;
        end

        % fi is custom tertinary operator
        cPs(((b-1)*T)+1:b*T) = fi(isHC, cP, 1-cP); % add controllability probabilities as a vector if is HC

        % Create state array
        sArr = repelem(cond, rep);
        % Suffle state array
        sArr = sArr(randperm(length(sArr)));
        % Create task environment
        env = TrialEnvironment(rP, cP, sArr, T);

        % Log sequences and vecs
        sSeq{b} = sArr;
        cVecs{b} = env.cVec;
        rVecs{b} = env.rVec;
        
        % Reset model parameters (except omega e.g.)
        m = m.resetModel();

        for t = 1:T
            tIdx = tIdx + 1;
            [s, cA, env] = env.presentTrial(); % get state, correct action and environment
            m = m.calcWeights(s); % calculate action weights

            m = m.calcProbs(s); % calc probabilities given new action weights

            % Create matrix with p(Go) values
            if isHC
                HCcell = updateProbGoMatrix(HCcell, hc, s, m);
            else
                LCcell = updateProbGoMatrix(LCcell, lc, s, m);
            end

            % Keep track of omega values
            if isO
                os(tIdx) = m.getOmega();
            end

            a = m.returnAction(s);
            aVec(b, t) = a;
            correct_aVec(b, t) = cA;
            r = env.getReward(t, s, a, isHC);
            rMat(b, t) = r;
            m = m.updateModel(r, s, a);
        end
    end
end

function matrix = updateProbGoMatrix(matrix, b, s, model)
    if isempty(matrix{b, s})
        matrix{b, s} = model.returnGoProb(s);
    else
        matrix{b, s}(end+1) = model.returnGoProb(s);
    end
end
