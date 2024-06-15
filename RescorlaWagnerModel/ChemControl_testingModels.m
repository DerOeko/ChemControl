% ChemControl_testingModels.m

% This is an interactive script---execute it step-by-step.
% It tests a series of models and averages their runs.
% Made for purposes of seeing whether the retrieved parameters return good
% estimates.
%
% EXPLANATION OF SETTINGS:
% dirs.root               = string, root directory of project.
%
% OUTPUTS:
% no outputs, just plots.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel
clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Set directories:

% Directories:
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.results = fullfile(dirs.root, 'Log', 'Behavior', 'Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
dirs.models         = fullfile(dirs.root, 'models');
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

% Input file
inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
data = load(inputFile).data;
fprintf('Loaded inputFile with %d subjects\n', size(data, 2));

nMod = length(dir(fullfile(dirs.models, "*.m")));

%% 01) LOAD AND SAVE:
%% 01a) Alternative A: Load model fitted with LAP:

parType = 'lap';

fname_mod = cell(nMod, 1);

for iMod = 1:nMod
    fname_mod{iMod} = fullfile(dirs.lap, sprintf('lap_mod%02d.mat', iMod));
end

fprintf("Load model %d fit with LAP\n", selMod);
fname       = load(fname_mod{selMod});
cbm         = fname.cbm;
subParam    = cbm.output.parameters;
nSub        = size(subParam, 1);
nParam      = size(subParam, 2);

% ----------------------------------------------------------------------- %
%% 01b) Alternative B: Load model fitted with HBI:

parType = 'hbi';
% extract
for iMod = 1:nMod
    % load
    fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d.mat', iMod));
    fprintf("Load model %d fit with HBI\n", iMod);
    load(fname_hbi)
    fprintf("Extract model %d fit with HBI\n", iMod);
    groupParams{iMod} = cbm.output.group_mean{:};
end

%% 02) SIMULATE
%% 02a) Run settings
nRuns = 500;
nTrials = 40;
nBlocks = 8;
nStates = 4;
HCmeans = zeros(nTrials/4, 4, nRuns);
fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

%% 02b) Run simulation and plot
for iMod = 1:nMod
    subplot(2, ceil(nMod/2), iMod);
    parameters = groupParams{iMod};
    for iRun = 1:nRuns
        subj = sim_subj();
        out = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", iMod));
        HCcell = out.HCcell;
        HCoccurrences = zeros(nTrials/4, 4);
        for s = 1:4
            for occurrence = 1:nTrials/4
                HCprobs = zeros(nBlocks/2, 1);
                for b = 1:nBlocks/2
                    HCprobs(b) = HCcell{b, s}(occurrence);
                end
                HCoccurrences(occurrence, s) = mean(HCprobs);
            end
        end
        HCmeans(:, :, iRun) = HCoccurrences; 
    end
    HCmeans = mean(HCmeans, 3);
    plotLearningCurves(HCmeans, sprintf("M%02d", iMod), true, fig1);
end


