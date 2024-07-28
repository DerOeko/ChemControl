function figHandle = plotStayAnalysis(data, figHandle)
nSub = size(data, 1);
nBlocks = size(data{1}.stimuli, 1);
nTrials = size(data{1}.stimuli, 2);
[hc_d, lc_d, yc_d] = extractControlTypeData(data);
dataTypes = {'hc_d', 'lc_d', 'yc_d'};
numTypes = length(dataTypes);

% Initialize P_stay for each control type
P_stay = cell(nSub, 4, numTypes); % Columns: [Go-Reward, Go-Punish, NoGo-Reward, NoGo-Punish]

% Loop through each control type
for dType = 1:numTypes
    data = eval(dataTypes{dType});
    % Initialize counters
    for i = 1:nSub
        P_stay{i, 1, dType} = [0, 0]; % [stay count, total count] for Go-Reward
        P_stay{i, 2, dType} = [0, 0]; % Go-Punish
        P_stay{i, 3, dType} = [0, 0]; % NoGo-Reward
        P_stay{i, 4, dType} = [0, 0]; % NoGo-Punish
    end

    % Process each subject
    for iSub = 1:nSub
        subj = data{iSub};
        actions = subj.actions;
        outcomes = subj.outcomes;
        stimuli = subj.stimuli;

        transformedOs = transformOutcomes(outcomes, stimuli);
        goActions = actions == 1;
        noGoActions = actions == 2;
        nBlocks = size(actions, 1);
        for iBlock = 1:nBlocks
            for iTrial = 1:nTrials-1 % Use nTrials-1 to avoid indexing out of bounds
                currentStimulus = stimuli(iBlock, iTrial);
                nextIdx = find(stimuli(iBlock, iTrial+1:end) == currentStimulus, 1, 'first') + iTrial;

                if ~isempty(nextIdx)
                    if goActions(iBlock, iTrial)
                        if transformedOs(iBlock, iTrial) % Go-Reward
                            P_stay{iSub, 1, dType}(2) = P_stay{iSub, 1, dType}(2) + 1;
                            if goActions(iBlock, nextIdx)
                                P_stay{iSub, 1, dType}(1) = P_stay{iSub, 1, dType}(1) + 1;
                            end
                        else % Go-Punish
                            P_stay{iSub, 2, dType}(2) = P_stay{iSub, 2, dType}(2) + 1;
                            if goActions(iBlock, nextIdx)
                                P_stay{iSub, 2, dType}(1) = P_stay{iSub, 2, dType}(1) + 1;
                            end
                        end
                    elseif noGoActions(iBlock, iTrial)
                        if transformedOs(iBlock, iTrial) % NoGo-Reward
                            P_stay{iSub, 3, dType}(2) = P_stay{iSub, 3, dType}(2) + 1;
                            if noGoActions(iBlock, nextIdx)
                                P_stay{iSub, 3, dType}(1) = P_stay{iSub, 3, dType}(1) + 1;
                            end
                        else % NoGo-Punish
                            P_stay{iSub, 4, dType}(2) = P_stay{iSub, 4, dType}(2) + 1;
                            if noGoActions(iBlock, nextIdx)
                                P_stay{iSub, 4, dType}(1) = P_stay{iSub, 4, dType}(1) + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

% Calculate the probabilities for each subject and control type
P_stay_prob = cell(nSub, 4, numTypes);
for dType = 1:numTypes
    for iSub = 1:nSub
        for j = 1:4
            if P_stay{iSub, j, dType}(2) > 0
                P_stay_prob{iSub, j, dType} = P_stay{iSub, j, dType}(1) / P_stay{iSub, j, dType}(2);
            else
                P_stay_prob{iSub, j, dType} = NaN; % Handle cases where there were no trials
            end
        end
    end
end

% Calculate the mean probabilities across subjects for each control type
mean_P_stay_prob = nan(numTypes, 4); % Rows: control types, Columns: [Go-Reward, Go-Punish, NoGo-Reward, NoGo-Punish]
sem_P_stay_prob = nan(numTypes, 4);  % Standard error of the mean

for dType = 1:numTypes
    for j = 1:4
        dataMatrix = cell2mat(P_stay_prob(:, j, dType));
        mean_P_stay_prob(dType, j) = nanmean(dataMatrix);
        sem_P_stay_prob(dType, j) = nanstd(dataMatrix) / sqrt(sum(~isnan(dataMatrix)));
    end
end

% Display mean results for each control type
for dType = 1:numTypes
    fprintf('%s:\n', dataTypes{dType});
    fprintf('Mean Go-Reward: %.2f\n', mean_P_stay_prob(dType, 1));
    fprintf('Mean Go-Punish: %.2f\n', mean_P_stay_prob(dType, 2));
    fprintf('Mean NoGo-Reward: %.2f\n', mean_P_stay_prob(dType, 3));
    fprintf('Mean NoGo-Punish: %.2f\n', mean_P_stay_prob(dType, 4));
end

% Plot the results as a bar plot with individual data points and error bars
figure;
hold on;
bar(mean_P_stay_prob', 'grouped');
xticklabels({'Go-Reward', '',  'Go-Punish', '', 'NoGo-Reward',  '',  'NoGo-Punish'});
yline(0.5, 'Color', '#080808', 'LineStyle', '--');

% Add individual subject data points
for dType = 1:numTypes
    for j = 1:4
        dataMatrix = cell2mat(P_stay_prob(:, j, dType));
        x = repmat(j + (dType-1)*0.22 - 0.22, size(dataMatrix)); % Adjust the position to align with bars
        scatter(x, dataMatrix, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.3);
    end
end

% Add error bars
for dType = 1:numTypes
    for j = 1:4
        errorbar(j + (dType-1)*0.22 - 0.22, mean_P_stay_prob(dType, j), sem_P_stay_prob(dType, j), 'k', 'LineStyle', 'none', 'LineWidth', 1);
    end
end

legend('high control', 'low control', 'yoked control', 'Location', 'Best');
xlabel('Condition');
ylabel('Mean Probability of Staying');
hold off;

end

