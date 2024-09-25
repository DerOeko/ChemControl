function [job] = ChemControl_cbm_sim(job)

%
% Perform model simulations for given model for given type of parameters.
%
% INPUTS:
% job 	      = structure with the following fields:
% .simType = string, type of simulation type, either 'realSim' (simulating
% based on real data) or 'modSim' (simulating based on synthetic data)
% .parType    = string, type of input parameters, either 'lap' (LaPlace approximation) or 'hbi' (Hierarchical Bayesian inference).
% .selMod       = integer, model number to be simulated.
%
% OUTPUTS:
% Save to disk.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.
% Inspired by code from Algermissen, J. et al. 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel
% clear all; close all; clc

% ----------------------------------------------------------------------- %

if ~isfield(job, 'simType')
    job.simType = 'modSim'; % modSim, realSim
end

if ~isfield(job, "parType")
    job.parType     = "lap";
end

if ~isfield(job, "selMod")
    selMod = 1;
    job.selMod      = selMod;
end

fprintf("Simulate for model %02d\n", job.selMod);
fprintf('Simulate data using method %s based on %s fits\n', ...
    job.simType, job.parType);
% ----------------------------------------------------------------------- %
%% Directories:

dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.results        = fullfile(dirs.root, 'Log', 'Behavior', 'Modelling_CBM');
if ~exist(dirs.results, "dir"); mkdir(dirs.results); end
dirs.lap            = fullfile(dirs.results, 'LAP_Results');
if ~exist(dirs.lap, "dir"); mkdir(dirs.lap); end
dirs.hbi            = fullfile(dirs.results, 'HBI_Results');
if ~exist(dirs.hbi, "dir"); mkdir(dirs.hbi); end
dirs.models         = fullfile(dirs.root, 'models');
if ~exist(dirs.models, "dir"); mkdir(dirs.models); end
dirs.sim            = fullfile(dirs.results, 'Simulations');
if ~exist(dirs.sim, "dir"); mkdir(dirs.sim); end

% Add paths:
addpath(fullfile(dirs.root, "modelSimulations"));

% ----------------------------------------------------------------------- %
%% Complete downstream settings:

job.nIter           = 100;

% Object name:
job.hbi_suffix      = sprintf("%02d", job.selMod);

% Complete name of HBI file:
job.hbi_name_mod    = fullfile(dirs.hbi, sprintf("hbi_mod_%s_allData.mat", job.hbi_suffix));

% Number of HBI model to inspect:
job.hbi_selMod = job.selMod;


% ----------------------------------------------------------------------- %
%% Define output file name:

outputFile = fullfile(dirs.sim,sprintf('%s_mod%02d_%s_iter%04d.mat', ...
    job.simType, job.selMod, job.parType, job.nIter));

% ----------------------------------------------------------------------- %
%% Start simulation:

if ~exist(outputFile, "file")

    % ------------------------------------------------------------------- %
    %% Load data:

    % Use original subject data --> original stimulus order
    fprintf("Load data\n");
    inputFile = fullfile(dirs.results, "ChemControl_cbm_inputData.mat");
    fdata = load(inputFile);
    data = fdata.data;
    job.nSub = length(data);

    % ------------------------------------------------------------------- %
    %% Load parameters:

    fprintf("Load parameters based on %s\n", job.parType)
    if strcmp(job.parType, 'lap')
        job.lap_name_mod = fullfile(dirs.lap, sprintf("lap_mod%02d_allData", job.selMod));
        fname = load(job.lap_name_mod);
        cbm = fname.cbm;
        allParam = cbm.output.parameters;
    elseif strcmp(job.parType, "hbi")
        fname = load(job.hbi_name_mod);
        cbm = fname.cbm;
        if length(cbm.output.parameters) == 1
            allParam = cbm.output.parameters{:};
        else
            allParam = cbm.output.parameters{job.hbi_selMod};
        end
    end

    % ------------------------------------------------------------------- %
    %% Simulate data:

    fprintf('Simulate %s based on %s parameters for model %02d with %d iterations\n', ...
        job.simType, job.parType, job.selMod, job.nIter)
    for iSub = 1:job.nSub

        % Extract subject data:
        fprintf("Start subject %03d\n", iSub)
        if strcmp(job.simType, 'realSim')
            subj = data{iSub};
        elseif strcmp(job.simType, 'modSim')
            subj = sim_subj(16, 80);
        else
            error("Unknown simulation type");
        end

        % Retrieve parameters:
        parameters = allParam(iSub, :);

        % Save per subject:
        sim.parameters{iSub}   = parameters; % add parameters
        sim.subj{iSub}         = subj; % add subject data

        % Iterate:
        for iIter = 1:job.nIter
            if job.nIter >= 100
                fraction_done = iIter/job.nIter;
                waitbar(fraction_done);
            end

            % Simulate:
            out = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", job.selMod));

            % Save results:
            sim.HCcell = out.HCcell;
            sim.LCcell = out.LCcell;
            sim.stimuli(iSub, iIter, :, :) = out.stimuli;
            sim.actions(iSub, iIter, :, :) = out.actions;
            sim.outcomes(iSub, iIter, :, :) = out.outcomes;
            sim.randHC(iSub, iIter, :, :) = out.randHC;
            sim.randLC(iSub, iIter, :, :) = out.randLC;

        end % End iIter
    end % End iSub

    % ------------------------------------------------------------------- %
    %% Save:

    fprintf('Save outputs\n')
    job.dirs = dirs;
    save(outputFile, 'sim', 'job', '-v7.3');

    fprintf('Done :------]\n');
% ----------------------------------------------------------------------- %
%% Otherwise load:
else
    
    fprintf("File %s already exists; \nload file\n", outputFile);
    load(outputFile);
    fprintf("File loaded :] \n");
end % End of if

end % END OF FUNCTION

    
