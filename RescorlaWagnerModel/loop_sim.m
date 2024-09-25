function loop_sim()

% loop_sim()
%
% Wrapper to call ChemControl_cbm_sim.m for model simulations
%
% INPUTS:
% None, set settings for simulations interactively.
%
% OUTPUTS:
% None, ChemControl_cbm_sim saves to disk.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.
% Inspired by code from Algermissen, J. et al. 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel
% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% Initialize root directory, add path:

dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
addpath(fullfile(dirs.root, "Log/Behavior/Modelling_CBM/Simulations"));

% ----------------------------------------------------------------------- %
%% Initialize job settings:

dirs.models         = fullfile(dirs.root, 'models');
nMod = length(dir(fullfile(dirs.models, "*.m")));
fprintf("Found %i models.\n", nMod)
simTypes = ["modSim"];
parTypes = ["lap" "hbi"];
job.simType = "modSim";
job.parType = "parType";

% ----------------------------------------------------------------------- %
%% Loop over models:

for i = 1:length(simTypes)
    for j = 1:length(parTypes)
        for iMod = 1:nMod
            job.simType = simTypes(i);
            job.parType = parTypes(j);
            job.selMod = iMod;
            fprintf('\n\n>>> Start simulating model %02d with simulation type %s and parameter type %s.\n', job.selMod, job.simType, job.parType);
            ChemControl_cbm_sim(job);
        end
    end
end

end % END OF FILE
