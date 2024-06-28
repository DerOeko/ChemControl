function figHandle = plotOmegas(averageOmegas6, averageOmegas7, cProbs, figHandle)
arguments
    averageOmegas6;
    averageOmegas7;
    cProbs;
    figHandle = [];
end

% Concatenate the averageOmegas cell arrays
averageOmegasList = [averageOmegas6, averageOmegas7];

% Initialize the figure
if isempty(figHandle)
    figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
else
    figure(figHandle);
end

% Iterate through each controllability schedule
nSchedules = size(averageOmegasList, 1);
for iSchedule = 1:nSchedules
    clf; % Clear the figure for the next plot

    % Get the average omegas for the current schedule for both models
    averageOmegas6 = averageOmegasList{iSchedule, 1};
    averageOmegas7 = averageOmegasList{iSchedule, 2};

    % Skip the schedule if either averageOmegas6 or averageOmegas7 is empty
    if isempty(averageOmegas6) || isempty(averageOmegas7)
        continue;
    end
    
    % Get the minimum and maximum values of averageOmegas
    minValue = min([averageOmegas6; averageOmegas7], [], 'all');
    maxValue = max([averageOmegas6; averageOmegas7], [], 'all');
    
    % Scale cProbs to the range of averageOmegas
    scaledOutcomeProbs = rescale(cProbs(iSchedule, :), minValue, maxValue);

    % Plot the data
    hold on
    plot(1:length(averageOmegas6), averageOmegas6, 'LineWidth', 2, 'DisplayName', 'Average Omegas 6')
    plot(1:length(averageOmegas7), averageOmegas7, 'LineWidth', 2, 'DisplayName', 'Average Omegas 7')
    plot(1:length(scaledOutcomeProbs), scaledOutcomeProbs, 'LineWidth', 3, 'LineStyle', '--', 'Color', '#AEAEAE', 'DisplayName', 'Scaled Outcome Probability')
    hold off

    xlabel('Trial number')
    ylabel('Omega value ([0, 1])')
    legend('Location', 'best')
    title(sprintf('Change of the Pavlovian Weighting \\omega over time for Schedule %d', iSchedule))

    % Wait for button press before proceeding to the next schedule
    disp('Press any key to continue to the next schedule...')
    waitforbuttonpress;
end

end
