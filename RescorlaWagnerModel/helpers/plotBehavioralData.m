function plotBehavioralData(out)
    % Plot the behavioral data from the model output

    % Unpack the output structure
    wsFreq = out.wsFreq;
    lsFreq = out.lsFreq;
    wsFrac = out.wsFrac;
    lsFrac = out.lsFrac;
    plotControl = out.plotControl;
    plotReward = out.plotReward;
    actions = out.actions;
    controllabilities = out.controllability;

    B = size(actions, 1);  % Number of blocks
    T = size(actions, 2);  % Number of trials per block

    figure;
    hold on;
    col = [0.1 0.7 0.1; 0.7 0.1 0.1];  % Colors for high and low control
    hs = [];  % Handles for strategy types
    hc = [];  % Handles for control types

    % Plotting Win-Stay and Lose-Shift frequencies by block with colors based on control
    for b = 1:B
        % Determine the color based on control level
        if controllabilities(b, 1) == 1
            idxColor = 1;  % High control (green)
        else
            idxColor = 2;  % Low control (red)
        end
        
        % Draw bars for Win-Stay and Lose-Shift
        h1 = bar(b-0.15, wsFreq(b), 0.3, 'FaceColor', col(idxColor,:));
        h2 = bar(b+0.15, lsFreq(b), 0.3, 'FaceColor', col(idxColor,:) + 0.3);  % Lighten for Lose-Shift
        
        % Collect handles for the legends
        if b == 1  % Only add once for the legend
            hs = [hs, h1];  % Win-Stay
            hs = [hs, h2];  % Lose-Shift
        end
        if idxColor == 1 && isempty(find(hc == 1, 1))
            hc(1) = h1;  % High control handle
        elseif idxColor == 2 && isempty(find(hc == 2, 1))
            hc(2) = h1;  % Low control handle
        end
    end
    
    hold off;
    title('Win-Stay and Lose-Shift Frequencies by Control Level');
    xlabel('Block');
    ylabel('Frequency');

    % Main legend for strategy types
    legend(hs, {'Win-Stay', 'Lose-Shift'}, 'Location', 'best');

    % Additional legend for control levels
    ah1 = axes('position', get(gca, 'position'), 'visible', 'off');
    legend(ah1, hc, {'High Control', 'Low Control'}, 'Location', 'bestoutside');
end
