function figHandle = plotShiftLoss(shiftMeans, iMod, figHandle)
arguments
    shiftMeans
    iMod = 1;
    figHandle = [];
end

% Initialize the figure
if isempty(figHandle)
    figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
else
    figure(figHandle); 
end

hold on
maxWins = size(shiftMeans, 2);
x = 1:maxWins;
barWidth = 0.2; % Width for each bar group

% Plot the results for each state with an offset
for iS = 1:4
    bar(x + (iS-2.5)*barWidth, shiftMeans(iS, :), barWidth, 'DisplayName', sprintf('State %d', iS));
end
hold off
xlabel('Number of Consecutive Wins');
ylabel('Weighted P(Shift | Lose after X Wins)');
title(sprintf("M%02d", iMod));
ylim([0 max(shiftMeans, [], 'all')]); % Set y-axis limits for better visualization
% Set xlim based on the last non-zero item in the first dimension of shiftMeans
lastNonZeroIndex = find(any(shiftMeans, 1), 1, 'last');
if isempty(lastNonZeroIndex)
    lastNonZeroIndex = maxWins; % Default to maxWins if all elements are zero
end
xlim([0 lastNonZeroIndex + 1]);
legend('Location','best')
end
