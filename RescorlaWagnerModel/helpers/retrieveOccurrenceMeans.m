function occurrenceMeans = retrieveOccurrenceMeans(probGoMatrix)
    T = length(probGoMatrix{1})*4;
    HC_B = size(probGoMatrix, 1);
    occurrenceMeans = NaN(T/4, 4);
    for s = 1:4
        for occurrence = 1:T/4
            occurrenceProbs = zeros(HC_B, 1);
            for b = 1:HC_B
                occurrenceProbs(b) = probGoMatrix{b, s}(occurrence);
            end
            occurrenceMeans(occurrence, s) = mean(occurrenceProbs);
        end
    end
end