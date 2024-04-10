function [averageHCprobGo, averageLCprobGo] = averageExperiment(numRuns, model, numTrialsInBlock, numBlocks, rewardProb, controllProb)
    averageHCprobGo = zeros(numTrialsInBlock/4, 4, numRuns);
    averageLCprobGo = zeros(numTrialsInBlock/4, 4, numRuns);

    for i = 1:numRuns
        [~, HCprobGoMatrix, LCprobGoMatrix, ~, ~, ~, ~, ~, ~] = runExperiment(model, numTrialsInBlock, numBlocks, rewardProb, controllProb);

        HCoccurrenceMeans = zeros(numTrialsInBlock/4, 4);
        for state = 1:4
            for occurrence = 1:numTrialsInBlock/4
                HCoccurrenceProbs = zeros(numBlocks/2, 1);
                for block = 1:numBlocks/2
                    HCoccurrenceProbs(block) = HCprobGoMatrix{block, state}(occurrence);
                end
                HCoccurrenceMeans(occurrence, state) = mean(HCoccurrenceProbs);
            end
        end

        LCoccurrenceMeans = zeros(numTrialsInBlock/4, 4);
        for state = 1:4
            for occurrence = 1:numTrialsInBlock/4
                LCoccurrenceProbs = zeros(numBlocks/2, 1);
                for block = 1:numBlocks/2
                    LCoccurrenceProbs(block) = LCprobGoMatrix{block, state}(occurrence);
                end
                LCoccurrenceMeans(occurrence, state) = mean(LCoccurrenceProbs);
            end
        end

        % Store the results for each run in the corresponding slices
        averageHCprobGo(:, :, i) = HCoccurrenceMeans;
        averageLCprobGo(:, :, i) = LCoccurrenceMeans;
    end

    % Calculate the mean across the third dimension (numRuns)
    averageHCprobGo = mean(averageHCprobGo, 3);
    averageLCprobGo = mean(averageLCprobGo, 3);
end