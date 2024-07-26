function figHandle = plotGoBias(data, label, figHandle)
    arguments
        data;
        label = "high control";
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
    
    % Calculate go bias
    P_Go_WinState = meanBlocks(1, :, :) - meanBlocks(3, :, :);
    
    % Mean across subjects, results in a 1 x Occurrences array
    meanGoBias = mean(P_Go_WinState, 3);
    
    % Standard error of the mean (SEM) across subjects
    semGoBias = std(P_Go_WinState, 0, 2) / sqrt(nSubs);

    % Plotting
    hold on;
    x = 1:nOccurrences;
    y = meanGoBias(:);
    sem = semGoBias(:);

    % Confidence bounds
    xconf = [x, x(end:-1:1)];
    yconf = [y + sem; flipud(y - sem)];

    % Plot confidence bounds
    fill(xconf, yconf, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Plot mean line
    plot(x, y, 'b', 'LineWidth', 2);

    xlabel('State Repetitions');
    xlim([1, nOccurrences]);
    ylabel('Go Bias');
    ylim([-1.0, 1.0]);
    yline(0, ":", 'LineWidth', 2, 'Color', '#AEAEAE');
    legend('', label , 'Location', 'best');

    grid on;
    hold off;
end
