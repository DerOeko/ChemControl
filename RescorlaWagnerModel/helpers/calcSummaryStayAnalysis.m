function stayAnalysisByControlType = calcSummaryStayAnalysis(out)
        % Store actions, outcomes and stimuli
        hc_actions = out.actions(find(out.controllability(:, 1) == 1), :);
        lc_actions = out.actions(find(out.controllability(:, 1) ~= 1), :);

        hc_outcomes = out.outcomes(find(out.controllability(:, 1) == 1), :);
        lc_outcomes = out.outcomes(find(out.controllability(:, 1) ~= 1), :);

        hc_stimuli = out.stimuli(find(out.controllability(:, 1) == 1), :);
        lc_stimuli = out.stimuli(find(out.controllability(:, 1) ~= 1), :);
    
        stayAnalysisByControlType = zeros(2, 4);

        nTrials = size(out.actions, 2);
        for i = 1:2
            if i == 1
                actions = hc_actions;
                outcomes = hc_outcomes;
                stimuli = hc_stimuli;
            elseif i == 2
                actions = lc_actions;
                outcomes = lc_outcomes;
                stimuli = lc_stimuli;
            end
            % Initialize stayAnalysis and counts for the current control type
            stayAnalysis = zeros(2, 2); % For Go/NoGo and Rewarded/Punished
            counts = zeros(2, 2); % For normalizing the results later

            for iBlock = 1:size(actions, 1)
                for iTrial = 1:nTrials-1 % Use nTrials-1 to avoid indexing out of bounds
                    currentStimulus = stimuli(iBlock, iTrial);
                    nextIdx = find(stimuli(iBlock, iTrial+1:end) == currentStimulus, 1, 'first') + iTrial;
                    if isempty(nextIdx)
                        continue; % Skip this iteration if no matching stimulus is found
                    end
                    currentOutcome = outcomes(iBlock, iTrial);
                    currentAction = actions(iBlock, iTrial);

                    nextAction = actions(iBlock, nextIdx);
                    
                    if currentOutcome == 1
                        if currentAction == 1
                            counts(1, 1) = counts(1, 1) + 1;
                            if nextAction == 1
                                stayAnalysis(1, 1) = stayAnalysis(1, 1) + 1;
                            end
                        else
                            counts(2, 1) = counts(2, 1) + 1;
                            if nextAction == 2
                                stayAnalysis(2, 1) = stayAnalysis(2, 1) + 1;
                            end
                        end
                    else % currentOutcome == -1
                        if currentAction == 1
                            counts(1, 2) = counts(1, 2) + 1;
                            if nextAction == 1
                                stayAnalysis(1, 2) = stayAnalysis(1, 2) + 1;
                            end
                        else
                            counts(2, 2) = counts(2, 2) + 1;
                            if nextAction == 2
                                stayAnalysis(2, 2) = stayAnalysis(2, 2) + 1;
                            end
                        end
                    end
                end
            end
            stayAnalysis = stayAnalysis ./ counts;	% Normalize by total counts
            stayAnalysis = reshape(stayAnalysis', 1, 4);
            stayAnalysisByControlType(i, :) = stayAnalysis(:);
        end
end

