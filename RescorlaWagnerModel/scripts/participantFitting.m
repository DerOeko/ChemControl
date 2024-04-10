%% Implementation of a Rescorla Wagner model Parameter retrieval
% Author: Samuel Nellessen, MCC
% samuelgerrit.nellessen@gmail.com/samuel.nellessen@ru.nl

% Clear previous sessions for a clean start
close all;
clear all;
clc;

%% Data Key:
% states: 1 - 8 (5-8 for low control, 1-4 for high control), int8
% controllArray: logical (1=high)
% actions: "None" or "space", logical (1 = Go)
% rewards: 10, 0, -10, int8
% rewardsBool: 1 = reward given/loss avoided in high control trials when correct response is given; in low control trials regardless response, 
% 0 = no reward/loss in high control trials when correct response is given; in low control trials regardless of response
% logical

%% Set up the Import Options and import the data
% Import the data
folderPath = "/project/3017083.01/behavioral study/data/raw/";
fileList = dir(fullfile(folderPath, "*.csv"));
participantData = struct();

for i = 1:numel(fileList)
    fileName = fileList(i).name;
    [~, name, ~] = fileparts(fileName);
    participantId = string(erase(name(1:7), "-"));

    filePath = fullfile(folderPath, fileName);
    opts = detectImportOptions(filePath); % Automatically detect options
    % Specify range and delimiter and variables
    opts.DataLines = [85, Inf];
    opts.Delimiter = ",";
    opts.MissingRule = "omitrow";
    opts.SelectedVariableNames = ["trialType", "miniBlock", "controllability", "randomReward", "feedback", "key_resp_keys"];
    data = readtable(filePath, opts);
    data = renamevars(data, ["trialType", "key_resp_keys"], ["state", "action"]);
    data.action = data.action == "space";
    data.state = int8(data.state);
    data.randomReward = data.randomReward == 1;
    data.controllability = data.controllability == "high";
    data.feedback = int8(data.feedback);
    data.miniBlock = int8(data.miniBlock);
    
    LCdata = data(~data.controllability, :);
    HCdata = data(data.controllability, :);
    participantData.(participantId).raw = data;
    participantData.(participantId).LCdata = LCdata;
    participantData.(participantId).HCdata = HCdata;
end

% Clear temporary variables
clear opts d

participants = fieldnames(participantData);
numStates = 4;

%% Other parameters, that stay constant and are not retrieved
numTrialsInBlock = 40;
numBlocks = 8;
rewardProb = 0.85;
controllProb = 0.8;
goBias = 0.3;
Qinit = zeros(4,2) + 0.5;
V = [0.5 -0.3 0.3 -0.3];
Vinit = [0.7 -0.3 0.3 -0.3];
pi = 0.3;

%% Get participant data
for index = 1:numel(participants)
    sub = string(participants(index));
    d = participantData.(sub).raw;
    controllabilityArray = zeros(numBlocks, 1);
    stimulusSequence = cell(numBlocks, 1);
    actions = zeros(numBlocks, numTrialsInBlock);
    correctActions = zeros(numBlocks, numTrialsInBlock);

    for block = 1:8
        miniBlockData = d(d.miniBlock == block, :);
        controllability = miniBlockData.controllability(1);
        controllabilityArray(block) = controllability;
        stimulusSequence{block} = miniBlockData.state;
        actions(block, :) = miniBlockData.action(:);
        correctActions(block, :) = miniBlockData.action(:);
end
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

function correctAction = getCorrectAction(state)
        % Get the correct action for a given state
        if mod(state, 4) == 1 || mod(state, 4) == 2
            correctAction = 1; % 'Go'
        else
            correctAction = 2; % 'NoGo'
        end
end
