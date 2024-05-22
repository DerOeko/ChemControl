function figHandle = plotParticipantCurves(participantData, isHC, figHandle)
    arguments
        participantData;
        isHC = true;
        figHandle = [];
    end
    if ~isHC
        % Define the function that will adjust LC state values
        adjustLCState = @(x) setfield(x, 'LCdata', ...
            setfield(x.LCdata, 'state', x.LCdata.state - 4));
        
        participantData = structfun(adjustLCState, participantData, 'UniformOutput', false);
    end


    controlString = fi(isHC, "High Control", "Low Control");
    run("config.m")
    
    if isempty(figHandle)
        figHandle = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    else
        figure(figHandle);
    end
    
    participants = fieldnames(participantData);
    numStates = 4;
    proportionGoResponses = zeros(numStates, 10, numel(participants));
    
    for i = 1:numel(participants)
        sub = string(participants(i));
        if isHC
            d = participantData.(sub).HCdata;
        elseif ~isHC
            d = participantData.(sub).LCdata;
        else
            disp("Specify HC flag");
        end

        for s = 1:numStates
                stateData = d(d.state == s, :);
                uniqueMiniblocks = unique(stateData.miniBlock);
                numMiniblocks = numel(uniqueMiniblocks);
                for m = 1:numMiniblocks
                    miniBlockData = stateData(stateData.miniBlock == uniqueMiniblocks(m, 1), :);
                    d.stateOccurrence(d.state == s & d.miniBlock == uniqueMiniblocks(m, 1)) = (1:10)';
                end
    
                for occ = 1:10
                    proportionGoResponses(s, occ, i) = size(d(d.state == s & d.action & d.stateOccurrence == occ, :), 1)/numMiniblocks;
                end
        end
        if isHC
            participantData.(sub).HCdata = d;
        elseif ~isHC
            participantData.(sub).LCdata = d;
        end
    end
    
    averageProportionGoResponses = mean(proportionGoResponses, 3);
    T = 40;
    
    hold on; % Allows multiple plots on the same figure
    for state = 1:4
        plot(1:T/4, averageProportionGoResponses(state, :)', 'LineWidth', 2); % Plotting mean probabilities for each state
    end
    xlabel('State Repetitions');
    xlim([1.0 T/4])
    ylabel('P(Go response | state)');
    ylim([0.0, 1.0])
    yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
    legend('Go to Win', 'Go to Avoid Loss', 'NoGo to Win', 'NoGo to Avoid Loss', 'Location', 'best');
    title_str = sprintf("Participant Data: \nMean P(Go|State) Across State Repetitions in\n%s Trials with N = %i", controlString, numel(participants));
    title(title_str);
    grid on
    
    hold off;
end