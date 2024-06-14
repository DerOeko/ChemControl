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

%% 00b) Settings:
nMod = length(dir(fullfile(dirs.models, '*.m')));
fprintf("Initialize settings with %d models\nSelected model %d\n", nMod, selMod);

selMod = 1;

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

% load
fname_hbi       = fullfile(dirs.hbi, sprintf('hbi_mod_%02d.mat', selMod));
fprintf("Load model %d fit with HBI\n", selMod);
load(fname_hbi)

% extract
fprintf("Extract model %d fit with HBI\n", selMod);
groupParam      = cbm.output.group_mean{:};
subParam        = cbm.output.parameters{:};
nSub            = size(subParam, 1);
nParam          = size(subParam, 2);

%% 02) SIMULATE

fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for iMod = 1:nMod
     subplot(3, 3, iMod);
     subj = sim_subj();
     parameters = [0.05 1.5 0.5 0.3 0.2 1.5 1.5];
     output = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", iMod));
     HCcell = output.HCcell;
     occurrenceMeans = retrieveOccurrenceMeans(HCcell);
     plotLearningCurves(occurrenceMeans, sprintf("ChemControl_mod%d_modSim", iMod), true, fig1);
end