function figHandle = plotOmegas(averageOmegasList, cProbs, figHandle)
arguments
    averageOmegasList;
    cProbs;
    figHandle = [];
end

% Initialize the figure
if isempty(figHandle)
    figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
else
    figure(figHandle);
end

% Iterate through each controllability schedule
nSchedules = size(averageOmegasList, 1);
nModels = size(averageOmegasList, 2); % Assuming each row has all averageOmegas

for iSchedule = 1:nSchedules
    clf; % Clear the figure for each schedule

    % Collect all the average omegas to determine the overall min and max
    allOmegas = [];
    for iModel = 1:nModels
        averageOmegas = averageOmegasList{iSchedule, iModel};
        if ~isempty(averageOmegas)
            allOmegas = [allOmegas; averageOmegas];
        end
    end

    % Skip the schedule if all averageOmegas are empty
    if isempty(allOmegas)
        continue;
    end

    % Prepare subplots for each model in a 2 by ceil(nModels/2) grid
    nCols = ceil(nModels / 4);
    disp(nCols)
    for iModel = 1:nModels
        subplot(4, nCols, iModel); % Create a subplot for each model

        averageOmegas = averageOmegasList{iSchedule, iModel};
        if ~isempty(averageOmegas)
            plot(1:length(averageOmegas), averageOmegas, 'LineWidth', 2, 'DisplayName', sprintf('Model %d', iModel + 5), 'Color', lines(1));
            ylabel('Omega value')
            hold on; % Keep plot for adding scaled probabilities

            % Customize each subplot with appropriate title
            title(sprintf('Model %d', iModel + 5));
            
            if iModel == nModels || iModel == nCols % Only add labels and legends on the last plots of each row
                xlabel('Trial number');
                legend('show', 'Location', 'best');
            end
        end
    end

    % Plot the outcome probabilities in all subplots
    for iModel = 1:nModels
        subplot(4, nCols, iModel);
        plot(1:length(cProbs(iSchedule, :)), cProbs(iSchedule, :), 'LineWidth', 3, 'LineStyle', '--', 'Color', '#AEAEAE', 'DisplayName', 'Scaled Outcome Probability')
    end

    % Add a super title for the entire schedule
    sgtitle(sprintf('Change of Pavlovian Weighting \\omega over time for Schedule %d', iSchedule));

    % Wait for button press before proceeding to the next schedule
    disp('Press any key to continue to the next schedule...')
    waitforbuttonpress;
end

end
