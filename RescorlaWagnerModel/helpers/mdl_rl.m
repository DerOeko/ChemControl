function x = mdl_rl(d, N, V)
%MDL_RL Computes Rissanen's Minimum Description Length for RL models.
% d = number of model parameters
% N = number of data points fitted (trials)
% V = negative log-likelihood of the model given the data

% Compute MDL
x = V * (1 + (d * log(N) / N));

end