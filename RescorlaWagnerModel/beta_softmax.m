function probabilities = beta_softmax(Q, beta)
    ex = exp(Q*beta);
    probabilities = ex/sum(ex);
end