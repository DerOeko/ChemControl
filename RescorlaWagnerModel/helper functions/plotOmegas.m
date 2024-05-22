function figHandle = plotOmegas(averageOmegasList, cProbs, cfg, figHandle)
arguments
    averageOmegasList;
    cProbs;
    cfg = Config();
    figHandle = [];
end
model_names = cfg.omegaModelNames;

if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
else
    figure(figHandle);
end
% Get the minimum and maximum values of averageOmegas
minValue = min(averageOmegasList, [], "all");
maxValue = max(averageOmegasList, [], "all");

% Scale cProbs to the range of averageOmegas
scaledOutcomeProbs = rescale(cProbs, minValue, maxValue);

% Plot the data
hold on
for o = 1:cfg.numOmegas
    plot(1:length(averageOmegasList(:, o)), averageOmegasList(:, o), "LineWidth", 2, "DisplayName", model_names{o})
end
plot(1:length(scaledOutcomeProbs), scaledOutcomeProbs, "LineWidth", 3, "LineStyle","--", 'Color', '#AEAEAE', "DisplayName", "Scaled Outcome Probability")
hold off
xlabel("Trial number")
ylabel("Omega value ([0, 1])")
legend('Location', 'best')
title("Change of the Pavlovian Weighting \omega over time")
end

