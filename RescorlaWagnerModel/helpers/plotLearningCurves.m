function figHandle = plotLearningCurves(occurrenceMeans, model_name, figHandle)
    arguments
        occurrenceMeans;
        model_name;
        figHandle = [];
    end
    sz = size(occurrenceMeans);
    T = sz(1) * sz(2);
    
    if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    else
        figure(figHandle);
    end
    
    hold on; % Allows multiple plots on the same figure
    
    for state = 1:4
        plot(1:T/4, occurrenceMeans(:, state)', 'LineWidth', 2); % Plotting mean probabilities for each state
    end

    xlabel('State Repetitions');
    xlim([1.0 T/4])
    ylabel('P(Go response | state)');
    ylim([0.0, 1.0])
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
    legend('G2W', 'G2A', 'NG2W', 'NG2A', 'Location', 'best');
    title(sprintf('%s', model_name));
    grid on
    hold off;
end