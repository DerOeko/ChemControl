clear all;
close all;
clc;

run('config.m');

%% Set up the Import Options and import the data
% Import the data
run("github_config.m")
fileList = dir(fullfile(folderPath, "*.csv"));
T = 40;
B = 8;
S = numel(fileList); 

data = cell(S, 1);
participantData = struct();
for s = 1:S
    fileName = fileList(s).name;
    [~, name, ~] = fileparts(fileName);
    participantId = string(erase(name(1:7), "-"));
    filePath = fullfile(folderPath, fileName);
    opts = detectImportOptions(filePath); % Automatically detect options
    opts.DataLines = [85, Inf];
    opts.Delimiter = ",";
    opts.MissingRule = "omitrow";
    opts.SelectedVariableNames = ["trialType", "miniBlock", "controllability", "randomReward", "feedback", "key_resp_keys"];
    
    d = readtable(filePath, opts);
    d = renamevars(d, ["trialType", "key_resp_keys"], ["state", "action"]);
    participantD = d;

    participantD.action = participantD.action == "space";
    participantD.state = int8(participantD.state);
    participantD.randomReward = participantD.randomReward == 1;
    participantD.controllability = participantD.controllability == "high";
    participantD.feedback = int8(participantD.feedback);
    participantD.miniBlock = int8(participantD.miniBlock);
    
    LCdata = participantD(~participantD.controllability, :);
    HCdata = participantD(participantD.controllability, :);
    participantData.(participantId).LCdata = LCdata;
    participantData.(participantId).HCdata = HCdata;

    d.controllability = d.controllability == "high";
    
    d.miniBlock = int8(d.miniBlock);

    % Convert actions from cell array of strings to a string array if necessary
    d.action = string(d.action);
    
    % Replace 'space' with '1' and 'None' with '2'
    d.action(d.action == "space") = "1";
    d.action(d.action == "None") = "2";
    
    % Convert the string array to integers
    d.action = str2double(d.action);  % This converts the string numbers to actual numeric (double) values
    
    % Process states
    d.state = double(d.state);
    d.state(d.state > 4) = d.state(d.state > 4) - 4;
    states = reshape(d.state, [T, B])';  % Reshape and then transpose
    
    % Process actions
    actions = reshape(d.action, [T, B])';  % Reshape and then transpose
    
    % Process outcomes
    d.feedback = double(d.feedback);
    outcomes = reshape(d.feedback, [T, B])';  % Reshape and then transpose
    
    % Store in a structured array
    d_struct = struct('states', states, 'actions', actions, 'outcomes', outcomes);
    data{s} = d_struct;
end

% Clear temporary variables
clear opts d

v = 6.25;

prior_model = struct('mean', zeros(2,1), 'variance', v);
fname_model = 'lap_model.mat';

prior_model_goBias = struct('mean', zeros(3, 1), 'variance', v);
fname_model_goBias = 'lap_model_goBias.mat';

prior_model_fixedPavlov = struct('mean', zeros(4, 1), 'variance', v);
fname_model_fixedPavlov = 'lap_model_fixedPavlov.mat';

prior_model_dynamicPavlov = struct('mean', zeros(4, 1), 'variance', v);
fname_model_dynamicPavlov = 'lap_model_dynamicPavlov.mat';

prior_model_fixedOmega = struct('mean', zeros(4, 1), 'variance', v);
fname_model_fixedOmega = 'lap_model_fixedOmega.mat';

prior_model_dynamicOmega1 = struct('mean', zeros(6, 1), 'variance', v);
fname_model_dynamicOmega1 = 'lap_model_dynamicOmega1.mat';

prior_model_dynamicOmega2 = struct('mean', zeros(7, 1), 'variance', v);
fname_model_dynamicOmega2 = 'lap_model_dynamicOmega2.mat';

cbm_lap(data, @model, prior_model, fname_model);
cbm_lap(data, @model_goBias, prior_model_goBias, fname_model_goBias);
cbm_lap(data, @model_fixedPavlov, prior_model_fixedPavlov, fname_model_fixedPavlov);
cbm_lap(data, @model_dynamicPavlov, prior_model_dynamicPavlov, fname_model_dynamicPavlov);
cbm_lap(data, @model_fixedOmega, prior_model_fixedOmega, fname_model_fixedOmega);
cbm_lap(data, @model_dynamicOmega1, prior_model_dynamicOmega1, fname_model_dynamicOmega1);
cbm_lap(data, @model_dynamicOmega2, prior_model_dynamicOmega2, fname_model_dynamicOmega2);

models = {@model, @model_goBias, @model_fixedPavlov, @model_dynamicPavlov, @model_fixedOmega, @model_dynamicOmega1, @model_dynamicOmega2};
fcbm_maps = {fname_model, fname_model_goBias, fname_model_fixedPavlov, fname_model_dynamicPavlov, fname_model_fixedOmega, fname_model_dynamicOmega1, fname_model_dynamicOmega2};
fname_hbi = 'hbi_model.mat';

cbm_hbi(data, models, fcbm_maps, fname_hbi);
cbm_hbi_null(data, fname_hbi)
d = load(fname_hbi);
cbm = d.cbm;
freq = cbm.output.model_frequency;
[~, k] = max(freq);
model_names = {'Model', 'GoBiasModel', 'FixedPavlovModel', 'DynamicPavlovModel', 'FixedOmegaModel', 'DynamicOmega1Model', 'DynamicOmega2Model'};

win_model = model_names{k};

switch win_model
    case 'Model'
        param_names = {'\epsilon', '\rho'};
        transform = {'sigmoid', 'exp'};

    case 'GoBiasModel'
        param_names = {'\epsilon', '\rho', 'goBias'};
        transform = {'sigmoid', 'exp', '@(x) x'};

    case 'FixedPavlovModel'
        param_names = {'\epsilon', '\rho', 'goBias', '\pi'};
        transform = {'sigmoid', 'exp', '@(x) x', '@(x) x'};

    case 'DynamicPavlovModel'
        param_names = {'\epsilon', '\rho', 'goBias', '\pi'};
        transform = {'sigmoid', 'exp', '@(x) x', '@(x) x'};

    case 'FixedOmegaModel'
        param_names = {'\epsilon', '\rho', 'goBias', '\omicron'};
        transform = {'sigmoid', 'exp', '@(x) x', 'sigmoid'};
    case 'DynamicOmega1Model'
        param_names = {'\epsilon', '\rho', 'goBias', '\omicron_{init}', '\alpha', '\kappa'};
        transform = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid', 'exp'};
    case "DynamicOmega2Model"
        param_names = {'\epsilon', '\rho', 'goBias', '\omega_{init}', '\alpha_{\Omega}','\beta_{\Omega}', '\thres_{\Omega}'};
        transform = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid', 'exp', '@tanh'};
    otherwise
        disp('Model not found.');
end

cbm_hbi_plot(fname_hbi, model_names, param_names, transform);
group_means = cbm.output.group_mean;
bestModelParameters = returnBestParameters(group_means{k}, transform);
fixedOmegaParameters = returnBestParameters(group_means{5}, {'sigmoid', 'exp', '@(x) x', 'sigmoid'});
dynamicQoppaParameters = returnBestParameters(group_means{6}, {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid', 'exp'});
dynamicOmegaParameters = returnBestParameters(group_means{7}, {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid', 'exp', '@tanh'});

% Display the best models parameters
disp(bestModelParameters);

bestModel = returnModel(win_model, bestModelParameters);
fixedOmegaModel = returnModel(model_names{5}, fixedOmegaParameters);
dynamicQoppaModel = returnModel(model_names{6}, dynamicQoppaParameters);
dynamicOmegaModel = returnModel(model_names{7}, dynamicOmegaParameters);

model_names = {'Generic Model', 'Go Bias Model', 'Fixed Motivational Bias Model', 'Dynamic Motivational Bias Model', "Fixed Omega Model", "Dynamic Omega Model with Qoppa", "Dynamic Omega Model with Omega"}; %"Dynamic Omega Model with Associability"};

best_model_name = model_names{k};
numRuns = numel(fieldnames(participantData));

fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
fig = figure(fig1);
averageOmegasList = zeros(numBlocks * T, 3);

% Plotting participants
subplot(2,5,1);
plotParticipantCurves(participantData, true, fig);
subplot(2,5,6);
plotParticipantCurves(participantData, false, fig);

% Plotting best model
[HCoccurrenceMeans, LCoccurrenceMeans, ~, outcomeProbs] = averageExperiment(numRuns, bestModel, T, numBlocks, rewardProb, controllProb, controllabilityArray);
subplot(2, 5, 2);
model_name = "Best fitting model: " + best_model_name;
plotLearningCurves(HCoccurrenceMeans, model_name, true, fig);
subplot(2, 5, 7);
plotLearningCurves(LCoccurrenceMeans, best_model_name, false, fig);

% Plotting fixed Omega Model
[HCoccurrenceMeans, LCoccurrenceMeans, averageOmegas, outcomeProbs] = averageExperiment(numRuns, fixedOmegaModel, T, numBlocks, rewardProb, controllProb, controllabilityArray);
averageOmegasList(:, 1) = averageOmegas;
subplot(2, 5, 3);
model_name = model_names{5};
plotLearningCurves(HCoccurrenceMeans, model_name, true, fig);
subplot(2, 5, 8);
plotLearningCurves(LCoccurrenceMeans, model_name, false, fig);

% Plotting dynamic Qoppa Model
[HCoccurrenceMeans, LCoccurrenceMeans, averageOmegas, outcomeProbs] = averageExperiment(numRuns, dynamicQoppaModel, T, numBlocks, rewardProb, controllProb, controllabilityArray);
averageOmegasList(:, 2) = averageOmegas;
subplot(2, 5, 4);
model_name = model_names{6};
plotLearningCurves(HCoccurrenceMeans, model_name, true, fig);
subplot(2, 5, 9);
plotLearningCurves(LCoccurrenceMeans, model_name, false, fig);

% Plotting dynamic Omega Model
[HCoccurrenceMeans, LCoccurrenceMeans, averageOmegas, outcomeProbs] = averageExperiment(numRuns, dynamicOmegaModel, T, numBlocks, rewardProb, controllProb, controllabilityArray);
averageOmegasList(:, 3) = averageOmegas;
subplot(2, 5, 5);
model_name = model_names{7};
plotLearningCurves(HCoccurrenceMeans, model_name, true, fig);
subplot(2, 5, 10);
plotLearningCurves(LCoccurrenceMeans, model_name, false, fig);

fig2 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
fig = figure(fig2);
plotOmegas(averageOmegasList, outcomeProbs, fig);
