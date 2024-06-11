function correctAction = getCorrectAction(state)
    % Get the correct action for a given state
    if state <= 2
        correctAction = 1; % 'Go'
    else
        correctAction = 2; % 'NoGo'
    end
end