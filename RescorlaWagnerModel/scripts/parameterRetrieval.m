%% Implementation of a Rescorla Wagner model Parameter retrieval
% Author: Samuel Nellessen, MCC
% samuelgerrit.nellessen@gmail.com/samuel.nellessen@ru.nl

% Clear previous sessions for a clean start
close all;
clear all;
clc;

% True hyperparameters
epsilon = 0.2;
beta = 3;
rho = 0.1;

%% Other parameters, that stay constant and are not retrieved
numTrialsInBlock = 160;
numBlocks = 28;
rewardProb = 0.85;
controllProb = 0.8;
goBias = 0.3;
Qinit = zeros(4,2) + 0.5;
V = [0.5 -0.3 0.3 -0.3];
Vinit = [0.7 -0.3 0.3 -0.3];
pi = 0.3;

%% Get simulated data
model = Model(epsilon, rho, beta, Qinit);
[blockInfo, HCprobGoMatrix, LCprobGoMatrix, stimulusSequence, outcomeVecs, successVecs, simulatedActions, correctActions, controllabilityArray] = runExperiment( ...
    model, ...
    numTrialsInBlock, ...
    numBlocks, ...
    rewardProb, ...
    controllProb);

%% Run parameter retrieval function
% Define parameter ranges
%paramRanges = createParamRanges();
paramRanges = {
    [epsilon - 0.05, epsilon, epsilon + 0.3, epsilon + 0.9], 
    [beta - 1, beta, beta + 5, beta + 10],
    [rho - 0.02, rho, rho + 0.1, rho + 0.5]
};
% Define model type
modelType = "Basic";
fprintf('Running parameter retrieval sanity check for %s Model...\n', modelType);
    fprintf('True parameters: epsilon = %.2f, beta = %.2f, rho = %.2f\n', epsilon, beta, rho);

optimalParams = retrieveParameters(simulatedActions, ...
    stimulusSequence, ...
    outcomeVecs, ...
    successVecs, ...
    correctActions, ...
    controllabilityArray, ...
    Qinit, ...
    Vinit, ...
    goBias, ...
    V, ...
    pi, ...
    modelType, ...
    paramRanges);
