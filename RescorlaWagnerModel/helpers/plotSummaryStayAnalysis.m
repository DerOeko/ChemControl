function figHandle = plotSummaryStayAnalysis(avgStayAnalysis, stderrStayAnalysis, model_name, figHandle)
    arguments
        avgStayAnalysis;
        stderrStayAnalysis;
        model_name;
        figHandle = [];
    end
    % Define labels for the categories
    categories = {'Go+', 'Go-', 'NoGo+', 'NoGo-', ''};
    controlTypes = {'HC', 'LC'};

    if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    else
        figure(figHandle);
    end
    
    dataToPlot = avgStayAnalysis';
    errorToPlot = stderrStayAnalysis';

    % Plot the bar graph
    barHandle = bar(dataToPlot);
    hold on;
    
    % Plot error bars
    numGroups = size(dataToPlot, 1); % Number of categories (4)
    numBars = size(dataToPlot, 2); % Number of control types (3)
    
    % Calculate the positions for error bars
    groupWidth = min(0.8, numBars/(numBars + 1.5));
    for i = 1:numBars
        x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
        errorbar(x, dataToPlot(:, i), errorToPlot(:, i), 'k', 'linestyle', 'none', 'CapSize', 5);
    end

    % Set x-axis labels
    set(gca, 'XTickLabel', categories, 'XTick', 1:numel(categories));



    % Add title and axis labels
    title(sprintf('%s', model_name));
    xlabel('Performed Action - Outcome');
    ylabel('P(stay)');

    % Set y-axis limits for better visibility (optional)
    ylim([0 1]);

    % Add grid for better readability (optional)
    grid on;

    % Adjust the color of the bars (optional)
    colormap('cool');
    yline(0.5, 'Color', '#080808', 'LineStyle', '--');
    % Add a legend to differentiate control types
    legend(controlTypes, 'Location', 'Best');
    hold off;
end

