function figHandle = plotBoxplots(HCoccurrenceMeans, LCoccurrenceMeans, model_name, figHandle)
    arguments
        HCoccurrenceMeans;
        LCoccurrenceMeans;
        model_name;
        figHandle = [];
    end
    
    run("config.m")
    
    if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    else
        figure(figHandle);
    end

    % Calculate Pavlovian bias
    GoToWin_GoToAvoid_HC = HCoccurrenceMeans(:, 1) - HCoccurrenceMeans(:, 2);
    NoGoToAvoid_NoGoToWin_HC = (HCoccurrenceMeans(:, 4)) - (HCoccurrenceMeans(:, 3));
    
    % Calculate Pavlovian bias
    GoToWin_GoToAvoid_LC = LCoccurrenceMeans(:, 1) - LCoccurrenceMeans(:, 2);
    NoGoToAvoid_NoGoToWin_LC = (LCoccurrenceMeans(:, 4)) - (LCoccurrenceMeans(:, 3));
    
    % Single matrix for boxplot
    biases = [GoToWin_GoToAvoid_HC, GoToWin_GoToAvoid_LC, NoGoToAvoid_NoGoToWin_HC, NoGoToAvoid_NoGoToWin_LC];
    boxplot(biases, 'Labels', {'GoToWin-GotoAvoid in HC', 'GoToWin-GotoAvoid in LC' 'NoGoToAvoid-NoGoToWin in HC', 'NoGoToAvoid-NoGoToWin in LC'});
    hold on;
    
    numDataPoints = size(biases, 1);
    % Create scatter plot for each group
    for i = 1:4
        scatter(repelem(i, numDataPoints), biases(:, i), 'd', 'filled');
    end
    
    hold off;
    grid on;
    title(sprintf('%s:\n Pavlovian Bias in High and Low Control Trials', model_name));
end