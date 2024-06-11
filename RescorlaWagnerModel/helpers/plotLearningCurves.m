function figHandle = plotLearningCurves(occurrenceMeans, model_name, isHC, figHandle)
    arguments
        occurrenceMeans;
        model_name;
        isHC = true;
        figHandle = [];
    end
    
    controlString = fi(isHC, "High Control", "Low Control");

    T = 40;
    
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
    title(sprintf('%s: \n P(Go|State) \n Across State Repetitions \nin %s Trials', model_name, controlString));
    grid on
    hold off;
end