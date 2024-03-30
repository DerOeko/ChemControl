function probabilities = betaSoftmax(actionWeights, beta)
    W_max = max(actionWeights);
    W_adj = actionWeights-W_max;
    ex = exp(W_adj*beta);
    probabilities = ex/sum(ex);
end