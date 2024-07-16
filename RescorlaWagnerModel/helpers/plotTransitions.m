function figHandle = plotTransitions(transitionData, figHandle)
    arguments
        transitionData;
        figHandle = [];
    end
    
    if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    else
        figure(figHandle);
    end
    transitions = fieldnames(transitionData);
    nTransitions = numel(transitions);
    nStates = 4;
    nOccurrences = 10;
    
    for i = 1:nTransitions
        subplot(3, 3, i);
        transition = transitions{i};
        data = transitionData.(transition);
        nBlocks = size(data, 2);
        transitionResponses = NaN(nStates, nOccurrences, nBlocks);
        for iBlock = 1:nBlocks
            blockData = data{iBlock};
            for iState = 1:nStates
                sIdx = find(blockData.stimuli(:) == iState);
                actions = blockData.actions(sIdx) == 1;
                transitionResponses(iState, 1:length(actions), iBlock) = actions;
            end
        end
        meanBlocks = squeeze(nanmean(transitionResponses, 3));
        % Calculate SEM
        semData = std(transitionResponses, 0, 3, 'omitnan') / sqrt(nBlocks);
        
        hold on;
        colors = lines(nStates);
        legendVec = ["G2W", "G2A", "NG2W", "NG2A"];
        for state = 1:nStates
            x = 1:nOccurrences;
            y = meanBlocks(state, :);
            sem = semData(state, :);
            
            % Confidence bounds
            xconf = [x, x(end:-1:1)];
            yconf = [y + sem, y(end:-1:1) - sem(end:-1:1)];
            
            % Plot confidence bounds
            fill(xconf, yconf, colors(state, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            
            % Plot mean line
            h1 = plot(x, y, 'LineWidth', 2, 'Color', colors(state, :), 'DisplayName', legendVec(state));
        end
        title(sprintf("%s with %i occurrences", transition, nBlocks));
        xlabel('Occurrences');
        ylabel('P(Go response | state)');
        ylim([0, 1]);
        xlim([1, nOccurrences]);
        legend('', 'G2W', '', 'G2A', '', 'NG2W', '', 'NG2A', 'Location', 'best');

        grid on;
        hold off;
    end
end