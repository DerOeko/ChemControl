function [HCmeans, LCmeans, averageOmegas, cPs] = averageExperiment(model, cfg)
    T = cfg.T;
    B = cfg.B;
    R = cfg.R;
    HCmeans = zeros(T/4, 4, R);
    LCmeans = zeros(T/4, 4, R);
    averageOmegas = zeros(T*B, R);


    for i = 1:R
        [HCcell, LCcell, ~, os, cPs] = runExperiment(model, cfg);

        HCoccurrences = zeros(T/4, 4);
        for s = 1:4
            for occurrence = 1:T/4
                HCprobs = zeros(B/2, 1);
                for b = 1:B/2
                    HCprobs(b) = HCcell{b, s}(occurrence);
                end
                HCoccurrences(occurrence, s) = mean(HCprobs);
            end
        end

        LCoccurrences = zeros(T/4, 4);
        for s = 1:4
            for occurrence = 1:T/4
                LCprobs = zeros(B/2, 1);
                for b = 1:B/2
                    LCprobs(b) = LCcell{b, s}(occurrence);
                end
                LCoccurrences(occurrence, s) = mean(LCprobs);
            end
        end
        
        averageOmegas(:, i) = os;
        % Store the results for each run in the corresponding slices
        HCmeans(:, :, i) = HCoccurrences;
        LCmeans(:, :, i) = LCoccurrences;

    end

    % Calculate the mean across the third dimension (numRuns)
    HCmeans = mean(HCmeans, 3);
    LCmeans = mean(LCmeans, 3);
    averageOmegas = mean(averageOmegas, 2);
end