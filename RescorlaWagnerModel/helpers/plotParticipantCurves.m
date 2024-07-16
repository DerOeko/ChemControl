function figHandle = plotParticipantCurves(data, figHandle)
    arguments
        data;
        figHandle = [];
    end
    
    if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    else
        figure(figHandle);
    end
    
    nStates = 4;
    nSubs = numel(data);
    nOccurrences = 10; % each state occurs 10 times per block
    nBlocks = size(data{1}.stimuli, 1);

    % Initialize the array to store "go" responses
    % Dimensions: States x Occurrences x Blocks x Subjects
    allResponses = NaN(nStates, nOccurrences, nBlocks, nSubs);
    for iSub = 1:nSubs
        d = data{iSub};
        for iBlock = 1:nBlocks
            for iState = 1:nStates
                sIdx = find(d.stimuli(iBlock, 1:40) == iState);
                actions = d.actions(iBlock, sIdx) == 1;
                allResponses(iState, 1:length(actions), iBlock, iSub) = actions;
            end
        end
    end

    % Mean across blocks, results in a States x Occurrences x Subjects array
    meanBlocks = squeeze(nanmean(allResponses, 3));
    
    % Mean across subjects, results in a States x Occurrences array
    meanSubjects = squeeze(nanmean(meanBlocks, 3));

    % Standard error of the mean (SEM) across subjects
    semSubjects = squeeze(nanstd(meanBlocks, [], 3) ./ sqrt(nSubs));

    hold on;
    colors = lines(nStates); % Get distinct colors for each state
    for state = 1:nStates
        x = 1:nOccurrences;
        y = meanSubjects(state, :);
        sem = semSubjects(state, :);

        % Confidence bounds
        xconf = [x, x(end:-1:1)];
        yconf = [y + sem, y(end:-1:1) - sem(end:-1:1)];
        
        % Plot confidence bounds
        p = fill(xconf, yconf, colors(state, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Plot mean line
        plot(x, y, 'LineWidth', 2, 'Color', colors(state, :));
    end
    
    xlabel('State Repetitions');
    xlim([1, nOccurrences]);
    ylabel('P(Go response | state)');
    ylim([0.0, 1.0]);
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE');
    legend('', 'G2W', '', 'G2A', '', 'NG2W', '', 'NG2A', 'Location', 'best');

    grid on;
    
    hold off;
end