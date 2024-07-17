function transformedOs = transformOutcomes(os, stimuli)

nBlocks = size(os, 1);
nTrials = size(os, 2);

for iBlock = 1:nBlocks
    for iTrial = 1:nTrials
        if os(iBlock, iTrial) == 1 && mod(stimuli(iBlock, iTrial), 2)
            transformedOs(iBlock, iTrial) = 1;
        elseif os(iBlock, iTrial) == 0 && ~mod(stimuli(iBlock, iTrial), 2)
            transformedOs(iBlock, iTrial) = 1;
        else
            transformedOs(iBlock, iTrial) = 0;
        end
    end
end

end

