clear all;
close all;
clc;

run('config.m');

epsilon_true = 0.3;
rho_true = 4;
goBias_true = 0.5;
pi_true = 0.5;
oi_true = 0.5;
o_true = 0.5;
alpha_true = 0.1;
kappa_true = 0.3;

Qinit = [0.5 -0.5 0.5 -0.5; 0.5 -0.5 0.5 -0.5]' * rho;
S = 40;   % Number of Subjects
B = numBlocks;
T = T;
data = cell(S, 1);

for s = 1:S
    modelClass = GoBiasModel(epsilon_true, rho_true, Qinit, goBias_true);
    [~, ~, ~, stateMat, ~, ~, actionMat, ~, ~, ~, ~, outcomeMat] = runExperiment(modelClass, ...
        T, ...
        numBlocks, ...
        rewardProb, ...
        controllProb, ...
        controllabilityArray);

    % States, Actions and Outcomes should be format (B, T)
    states = reshape(cell2mat(stateMat), [B, T]);
    actions = reshape(actionMat, [B, T]);
    outcomes = reshape(outcomeMat, [B, T]);
    d_struct = struct('states', states, 'actions', actions, 'outcomes', outcomes);
    data{s} = d_struct;

end


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

prior_model_dynamicOmega = struct('mean', zeros(6, 1), 'variance', v);
fname_model_dynamicOmega = 'lap_model_dynamicOmega.mat';

cbm_lap(data, @model, prior_model, fname_model);
cbm_lap(data, @model_goBias, prior_model_goBias, fname_model_goBias);
cbm_lap(data, @model_fixedPavlov, prior_model_fixedPavlov, fname_model_fixedPavlov);
cbm_lap(data, @model_dynamicPavlov, prior_model_dynamicPavlov, fname_model_dynamicPavlov);
cbm_lap(data, @model_fixedOmega, prior_model_fixedOmega, fname_model_fixedOmega);
cbm_lap(data, @model_dynamicOmega, prior_model_dynamicOmega, fname_model_dynamicOmega);

models = {@model, @model_goBias, @model_fixedPavlov, @model_dynamicPavlov, @model_fixedOmega, @model_dynamicOmega};
fcbm_maps = {fname_model, fname_model_goBias, fname_model_fixedPavlov, fname_model_dynamicPavlov, fname_model_fixedOmega, fname_model_dynamicOmega};
fname_hbi = 'hbi_model.mat';

cbm_hbi(data, models, fcbm_maps, fname_hbi);
cbm_hbi_null(data, fname_hbi)
d = load(fname_hbi);
cbm = d.cbm;
freq = cbm.output.model_frequency;
[~, k] = max(freq);
model_names = {'Model', 'GoBiasModel', 'FixedPavlovModel', 'DynamicPavlovModel', 'FixedOmegaModel', 'DynamicOmegaModel'};

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

    case 'DynamicOmegaModel'
        param_names = {'\epsilon', '\rho', 'goBias', '\omicron_{init}', '\alpha', '\kappa'};
        transform = {'sigmoid', 'exp', '@(x) x', 'sigmoid', 'sigmoid', 'exp'};
    otherwise
        disp('Model not found.');
end

cbm_hbi_plot(fname_hbi, model_names, param_names, transform);


