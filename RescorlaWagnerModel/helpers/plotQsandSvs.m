function figHandle = plotQsandSvs(averageQsCell, averageSvsCell, cProbs, modelsWithQs, figHandle)
arguments
    averageQsCell;
    averageSvsCell;
    cProbs;
    modelsWithQs = [];
    figHandle = [];
end

% Initialize the figure
if isempty(figHandle)
    figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
else
    figure(figHandle);
end

% Iterate through each controllability schedule
nModels = size(averageQsCell, 2); % Assuming each row has all averageOmegas

% Prepare subplots for each model in a 2 by ceil(nModels/2) grid
nCols = ceil(nModels / 4);
for iModel = 1:nModels
    modelName = modelsWithQs(iModel);
    subplot(4, nCols, iModel); % Create a subplot for each model
    
    averageQs = averageQsCell{iModel};
    averageSvs = averageSvsCell{iModel};
    if ~isempty(averageQs)
        hold on; % Keep plot for adding scaled probabilities
        plot(1:length(averageQs), averageQs, 'LineWidth', 2, 'DisplayName', "Q value");
        plot(1:length(averageSvs), averageSvs, 'LineWidth', 2, 'DisplayName', "V value");

        ylabel('Predicted value')
        
    
        % Customize each subplot with appropriate title
        title(sprintf('%s: %.2f', modelName, corr(averageQs, averageSvs)));
        
        if iModel == nModels || iModel == nCols % Only add labels and legends on the last plots of each row
            xlabel('Trial number');
            legend('show', 'Location', 'best');
        end
    end
end

% Plot the outcome probabilities in all subplots
for iModel = 1:nModels
subplot(4, nCols, iModel);
plot(1:length(cProbs(1, :)), cProbs(1, :), 'LineWidth', 3, 'LineStyle', '--', 'Color', '#AEAEAE', 'DisplayName', 'Outcome Probability')
end

% Add a super title for the entire schedule
sgtitle('Change of value estimation over time for first schedule and first stimulus');


end
