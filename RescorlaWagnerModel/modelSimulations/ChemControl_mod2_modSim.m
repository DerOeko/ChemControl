function [out] = ChemControl_mod2_modSim(parameters, subj)
% Standard Q-learning model with delta learning rule.
    
    % ----------------------------------------------------------------------- %
    %% Retrieve parameters:
    ep = sigmoid(parameters(1));
    rho = exp(parameters(2));
    gB = parameters(3);
    % ----------------------------------------------------------------------- %

    %% Unpack data:

    % Extract task features:
    stimuli = subj.stimuli; 
    controllabilities = subj.controllability; % 1 or 0
    randomRewards = subj.randomReward; % 1 or 0
    randHCs = subj.randHC; % 1, 0, 2
    randLCs = subj.randLC; % 1, 0, 2
    cali_stimuli = subj.cali_stimuli;
    cali_randHC = subj.cali_randHC;
    cali_randLC = subj.cali_randLC;
    cali_randRewards = subj.cali_randReward;

    % Data dimensions:
    B = size(stimuli, 1); % Number of blocks
    T = size(stimuli, 2); % Number of trials
    S = 4; % Number of stimuli

    % Store outputs
    HCcell = cell(B/2, S); % Store go probs in HC
    LCcell = cell(B/2, S); % Store go probs in LC
    actions = zeros(B, T);
    outcomes = zeros(B, T);

    q0 = [0.5 -0.5 0.5 -0.5];
    hc = 0;
    lc = 0;

    %% Run calibration block
    rewardLossCounter = zeros([1, 2]);
    q_g = q0 * rho;
    q_ng = q0 * rho;
    w_g = q0 * rho;
    w_ng = q0 * rho;
    isHC = 1;
    for t = 1:T
        s = cali_stimuli(t);
        randHC = cali_randHC(t);
        randLC = cali_randLC(t);
        isRewarded = cali_randRewards(t);

        w_g(s) = q_g(s) + gB;
        w_ng(s) = q_ng(s);
        p1 = stableSoftmax(w_g(s), w_ng(s));

        a = returnAction(p1);
        o = returnReward(s, a, isHC, randLC, randHC, isRewarded);

        if a==1
            q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
        elseif a==2
            q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
        end

        counter = updateRewardLossCounter(s, o);

        rewardLossCounter = rewardLossCounter + counter;
    end
   
    M = rewardLossCounter/(T/2);
    averageRewardRate = M(1);
    averageLossRate = M(2);
    numRewarded = round(averageRewardRate * T); % number of rewarded trials
    numAvoided = round(averageLossRate * T);

    for b = 1:B
        switch controllabilities(b, 1)
            case 1
                isHC = true;
                isLC = false;
                isYoked = false;
            case 0
                isHC = false;
                isLC = true;
                isYoked = false;
            case 2
                isHC = false;
                isLC = true;
                isYoked = true;

                rewardedVec = [ones(1, numRewarded) zeros(1, T-numRewarded, 1)];
                rewardedVec = rewardedVec(randperm(length(rewardedVec)));
                avoidedVec = [ones(1, numAvoided) zeros(1, T-numAvoided, 1)];
                avoidedVec = avoidedVec(randperm(length(avoidedVec)));
        end
        hc = hc + isHC;
        lc = lc + isLC;
        
        q_g = q0 * rho;
        q_ng = q0 * rho;
        w_g = q0 * rho;
        w_ng = q0 * rho;

        for t = 1:T
            s = stimuli(b, t);
            isWinState = mod(s, 2);
            randHC = randHCs(b, t); % outcome matters (1, 0, 2)
            randLC = randLCs(b, t); % outcome doesn't matter (1, 0, 2)
            if ~isYoked
                isRewarded = randomRewards(b, t);
            elseif isWinState && isYoked
                isRewarded = rewardedVec(t);
            elseif ~isWinState && isYoked
                isRewarded = avoidedVec(t);
            end

            w_g(s) = q_g(s) + gB;
            w_ng(s) = q_ng(s);

            p1 = stableSoftmax(w_g(s), w_ng(s));
            
            if isHC
                HCcell{hc, s}(end+1) = p1;
            else
                LCcell{lc, s}(end+1) = p1;
            end
            
            a = returnAction(p1);
            o = returnReward(s, a, isHC, randLC, randHC, isRewarded);
            actions(b, t) = a;
            outcomes(b, t) = o;
            if a==1
                q_g(s) = q_g(s) + ep * (rho * o - q_g(s));
            elseif a==2
                q_ng(s) = q_ng(s) + ep * (rho * o - q_ng(s));
            end
        end
    end
    % ----------------------------------------------------------------------- %
    %% Win stay-lose shift:

    win_stay_count = zeros(B, T);
    lose_shift_count = zeros(B, T);
    total_wins = zeros(B, 1);
    total_losses = zeros(B, 1);
    
    for b = 1:B
        for t = 1:T
            s = stimuli(b, t);
            current_action = actions(b, t);
            current_outcome = outcomes(b, t);
    
            % Find the next index after 't' where the stimulus 's' appears
            next_idx = find(stimuli(b, t+1:end) == s, 1) + t;
    
            if ~isempty(next_idx) % Ensure there is a next occurrence
                next_action = actions(b, next_idx);
    
                if current_outcome == 1 % Win condition
                    total_wins(b) = total_wins(b) + 1;
                    if current_action == next_action
                        win_stay_count(b, t) = 1;
                    end
                elseif current_outcome == -1 % Lose condition
                    total_losses(b) = total_losses(b) + 1;
                    if current_action ~= next_action
                        lose_shift_count(b, t) = 1;
                    end
                end
            end
        end
    end
    
    % Calculate frequencies for each block and overall
    wsFreq = sum(win_stay_count, 2) ./ total_wins;
    lsFreq = sum(lose_shift_count, 2) ./ total_losses;
    
    % Aggregate frequencies across blocks
    owsFreq = sum(sum(win_stay_count)) / sum(total_wins);
    olsFreq = sum(sum(lose_shift_count)) / sum(total_losses);

    % Cumulative win-stay-lose-shift analysis
    cumulative_win_stays = zeros(B, T);
    cumulative_lose_shifts = zeros(B, T);
    cumulative_wins = zeros(B, T);
    cumulative_losses = zeros(B, T);
    
    for b = 1:B
        for t = 1:T
            if t > 1
                cumulative_win_stays(b, t) = cumulative_win_stays(b, t-1) + win_stay_count(b, t);
                cumulative_lose_shifts(b, t) = cumulative_lose_shifts(b, t-1) + lose_shift_count(b, t);
                cumulative_wins(b, t) = cumulative_wins(b, t-1) + (outcomes(b, t) == 1);
                cumulative_losses(b, t) = cumulative_losses(b, t-1) + (outcomes(b, t) == -1);
            else
                cumulative_win_stays(b, t) = win_stay_count(b, t);
                cumulative_lose_shifts(b, t) = lose_shift_count(b, t);
                cumulative_wins(b, t) = (outcomes(b, t) == 1);
                cumulative_losses(b, t) = (outcomes(b, t) == -1);
            end
        end
    end
    
    wsFrac = cumulative_win_stays ./ cumulative_wins;
    lsFrac = cumulative_lose_shifts ./ cumulative_losses;
    
    wsFreqHC = zeros(4, 1);
    wsFreqLC = zeros(4, 1);
    lsFreqHC = zeros(4, 1);
    lsFreqLC = zeros(4, 1);
    total_winsHC = zeros(4, 1);
    total_lossesHC = zeros(4, 1);
    total_winsLC = zeros(4, 1);
    total_lossesLC = zeros(4, 1);

    for b = 1:B
        if controllabilities(b, 1) == 1
            wsFreqHC(b) = wsFreq(b);
            lsFreqHC(b) = lsFreq(b);
            total_winsHC = total_wins(b);
            total_lossesHC = total_losses(b);
        else
            wsFreqLC(b) = wsFreq(b);
            lsFreqLC(b) = lsFreq(b);
            total_winsLC = total_wins(b);
            total_lossesLC = total_losses(b);
        end    
    end
    
    % Initialize vectors
    maxWins = 40; % Maximum possible wins
    %shiftAfterLossCounts = zeros(B, maxWins);
    shiftAfterLossCounts = zeros(4, maxWins);
    totalConsecutiveWinsCounts = zeros(4, maxWins);
    weightedProbShiftAfterLoss = zeros(4, maxWins);
    S = 4; % Number of stimuli

    for iS = 1:S
        for b = 1:B
            consecutiveWins = 0;
            for t = 1:T - 1
                if stimuli(b, t) ~= iS
                    continue
                end
                s = stimuli(b, t);
                o = outcomes(b, t);
                isWinState = mod(s, 2);
                if (isWinState && o == 1) || (~isWinState && o == 0)
                    consecutiveWins = consecutiveWins + 1;
                else
                    if consecutiveWins > 0
                        % Count the total number of consecutive wins
                        totalConsecutiveWinsCounts(iS, consecutiveWins) = totalConsecutiveWinsCounts(iS, consecutiveWins) + 1;                        
                        nextStateTrial = find(stimuli(b, t+1:end) == iS, 1, 'first') +t;
                        if ~isempty(nextStateTrial) && actions(b, t) ~= actions(b, nextStateTrial)
                            shiftAfterLossCounts(iS, consecutiveWins) = shiftAfterLossCounts(iS, consecutiveWins) + 1;
                        end
                    end
                    consecutiveWins = 0;
                end
            end
        end
    end
   
    % Calculate probabilities
    probShiftAfterLoss = shiftAfterLossCounts ./ totalConsecutiveWinsCounts;
    
    % Handle division by zero (NaN values)
    probShiftAfterLoss(isnan(probShiftAfterLoss)) = 0;

    % Weight the probabilities by their occurrences
    weightedProbShiftAfterLoss = probShiftAfterLoss .* totalConsecutiveWinsCounts ./ sum(totalConsecutiveWinsCounts);
    
    % Controllability array (for plotting)
    plotControl = zeros(B, T);
    for b = 1:B
        if controllabilities(b, 1) == 1
            plotControl(b, :) = 0.8;
        else
            plotControl(b, :) = 0.2;
        end
    end

    % Reward array (for plotting)
    plotReward = zeros(B, T);
    for b = 1:B
        if controllabilities(b, 1) == 1 || controllabilities(b, 1) == 0
            plotReward(b, :) = 0.8;
        else
            plotReward(b, :) = (averageLossRate + averageRewardRate) /2;
        end
    end

    % ----------------------------------------------------------------------- %
    %% Save as output object:
    out.HCcell = HCcell;
    out.LCcell = LCcell;
    out.randHC = randHCs;
    out.randLC = randLCs;
    out.stimuli = stimuli;
    out.controllability = controllabilities;
    out.randomRewards = randomRewards;
    out.actions = actions;
    out.outcomes = outcomes;
    out.wsFreq = wsFreq;
    out.lsFreq = lsFreq;
    out.owsFreq = owsFreq;
    out.olsFreq = olsFreq;
    out.wsFrac = wsFrac;
    out.lsFrac = lsFrac;
    out.plotControl = plotControl;
    out.plotReward = plotReward;
    out.arr = averageRewardRate;
    out.alr = averageLossRate;
    out.probShiftAfterLoss = probShiftAfterLoss;
    out.weightedProbShiftAfterLoss = weightedProbShiftAfterLoss;

end 
