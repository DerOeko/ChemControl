function o = returnReward(s, a, isHC, randLC, randHC, isRewarded)
    isWinState = mod(s, 2);

    isCorrectAction = a == getCorrectAction(s);

    if isHC % high control block
            if randHC == 1 % outcome matters
                if isCorrectAction % correct action
                    if isRewarded % correct actions are rewarded
                        if isWinState % win state
                            o = 1;
                        else 
                            o = 0; % loss state (neutral outcome, if correct action)
                        end % End of win state check
                    else % no reward even if correct response
                        if isWinState
                            o = 0;
                        else
                            o = -1;
                        end
                    end
                else % if not correct response
                    if isRewarded
                        if isWinState
                            o = 0;
                        else
                            o = -1;
                        end
                    else
                        if isWinState
                            o = 1;
                        else
                            o = 0;
                        end
                    end
                end % end of isCorrectAction
            elseif randHC == 0 % if outcome doesn't matter
                if isRewarded
                    if isWinState
                        o = 1;
                    else
                        o = 0;
                    end
                else
                    if isWinState
                        o = 0;
                    else
                        o = -1;
                    end
                end
            end % end of isNoControl
        else % isLowControl
            if randLC == 0
                if isCorrectAction % correct action
                    if isRewarded % correct actions are rewarded
                        if isWinState % win state
                            o = 1;
                        else
                            o = 0; % loss state (neutral outcome, if correct action)
                        end % End of win state check
                    else % no reward even if correct response
                        if isWinState
                            o = 0;
                        else
                            o = -1;
                        end
                    end
                else % if not correct response
                    if isRewarded
                        if isWinState
                            o = 0;
                        else
                            o = -1;
                        end
                    else
                        if isWinState
                            o = 1;
                        else
                            o = 0;
                        end
                    end
                end % end of isCorrectAction
            elseif randLC == 1
                if isRewarded
                    if isWinState
                        o = 1;
                    else
                        o = 0;
                    end
                else
                    if isWinState
                        o = 0;
                    else
                        o = -1;
                    end
                end % end of isRewarded
            end % end of randLC
        end % end of isHC
end

