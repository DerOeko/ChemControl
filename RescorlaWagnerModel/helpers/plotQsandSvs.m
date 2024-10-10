function figHandle = plotQsandSvs(s, averageQs, averageSvs, cProbs, modelsWithQs, figHandle)
arguments
    s;
    averageQs;
    averageSvs;
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
nModels = size(averageOmegasCell, 2); % Assuming each row has all averageOmegas

% Prepare subplots for each model in a 2 by ceil(nModels/2) grid
nCols = ceil(nModels / 4);
for iModel = 1:nModels
    modelName = modelsWithOmega(iModel);
    subplot(4, nCols, iModel); % Create a subplot for each model
    
    averageOmegas = averageOmegasCell{iModel};
    if ~isempty(averageOmegas)
        plot(1:length(averageOmegas), averageOmegas, 'LineWidth', 2, 'DisplayName', modelName, 'Color', lines(1));
        ylabel('Omega value')
        hold on; % Keep plot for adding scaled probabilities
    
        % Customize each subplot with appropriate title
        title(modelName);
        
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
sgtitle('Change of Pavlovian Weighting \omega over time for  first schedule');


end
