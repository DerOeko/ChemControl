function counter = updateRewardLossCounter(s, o)
    isWinState = mod(s, 2);

    if isWinState && o == 1
        counter = [1 0];
    elseif ~isWinState && o == 0
        counter = [0 1];
    else
        counter = [0 0];
    end
end

