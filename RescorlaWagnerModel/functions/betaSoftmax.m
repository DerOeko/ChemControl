function p = betaSoftmax(w, beta)
    p = exp(w*beta)/sum(exp(w*beta));
end