function probabilities = beta_softmax(Q, beta)
    Q_max = max(Q);
    Q_adj = Q-Q_max;
    ex = exp(Q_adj*beta);
    probabilities = ex/sum(ex);
end