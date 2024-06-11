% ChemControl_cbm_prepareData.m
% 
% Execute this script to prepare the behavioral data for an appropriate format for the CBM toolbox.
% Mind adjusting the root directory.
%
% INPUTS:
% none.
%
% OUTPUTS:
% Saves all data as ChemControl_cbm_inputData.mat to indicated directory.
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.
% Inspired by code from Algermissen, J. et al. 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
warning('off','all')

%% Initialize directories:

run("github_config.m") % Contains path to actual data folder; adjust with your onw path
dirs = [];
dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.behav = fullfile(folderPath);
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
if ~exist(dirs.target, 'dir'); mkdir(dirs.target); end

% ----------------------------------------------------------------------- %
%% Settings:

fprintf("Initialize settings\n")

% Fixed input settings:
nCueConditions = 4; % number of cue conditions, here 4 for Go2Win, Go2Avoid, NoGo2Win, NoGo2Avoid
nExemplar = 10; % how often each cue is shown 
nBlocks = 8; % how many blocks in one experiment
nTrials = 40; % how many trials in a block

% Downstream settings:

% ----------------------------------------------------------------------- %
%% Retrieve subject IDs:

subFind = dir(fullfile(dirs.behav, "*.csv"));
subFind = {subFind.name};
subList = nan(numel(subFind), 1);

for iSub = 1:numel(subFind)
    subList(iSub) = str2double(string(extractBetween(subFind(iSub), "sub-", "_Robo")));
end

sID = sort(subList);
nSub = numel(sID);

% ----------------------------------------------------------------------- %
%% Load rawData:

% Preallocate a cell array to store all the data
rawData = cell(nSub, 1);

% Load raw data
fprintf("Load raw data\n");

for iSub = 1:nSub
    fprintf('Load data for subject %03d\n', iSub);
    
    % Define file path
    filePath = fullfile(folderPath, subFind{iSub});
    
    % Define import options
    opts = detectImportOptions(filePath);
    opts.DataLines = [85, Inf]; % Start reading from line 85
    opts.Delimiter = ","; % CSV delimiter
    opts.MissingRule = "omitrow"; % Omit rows with missing data
    opts.SelectedVariableNames = ["trialType", "miniBlock", "controllability", "randHC", "randLC", "randomReward", "feedback", "key_resp_keys"];
    
    % Try to load the data
    try
        data = readtable(filePath, opts);
    catch ME
        fprintf('Error loading data for subject %03d: %s\n', iSub, ME.message);
        continue;
    end
    
    % Rename variables for consistency
    data = renamevars(data, ["trialType", "key_resp_keys", "feedback"], ["stimuli", "actions", "outcomes"]);
    
    % Convert data types and values
    data.stimuli(data.stimuli > 4, :) = data.stimuli(data.stimuli > 4) - 4; 
    data.miniBlock = int8(data.miniBlock);
    data.controllability = strcmp(data.controllability, 'high');
    data.randHC = int8(data.randHC);
    data.randLC = int8(data.randLC);
    data.randomReward = data.randomReward == 1;
    % Convert actions from cell array of strings to a string array if necessary
    data.actions = string(data.actions);
    
    % Replace 'space' with '1' and 'None' with '2'
    data.actions(data.actions == "space") = "1";
    data.actions(data.actions == "None") = "2";
    % Convert the string array to integers
    data.actions = str2double(data.actions);  % This converts the string numbers to actual numeric (double) values
    data.outcomes = double((data.outcomes == 10) - (data.outcomes == -10)); % Convert outcomes to -1, 0, 1
    
    % Store the data for this subject in the cell array
    rawData{iSub} = data;
end

% ----------------------------------------------------------------------- %
%% Reshape into structure:

fprintf('Reshape data\n');

% Initialize empty structure:
data = struct([]);

for iSub = 1:nSub
    data{iSub}.stimuli = reshape(rawData{iSub}.stimuli, [nTrials, nBlocks])';
    data{iSub}.actions = reshape(rawData{iSub}.actions, [nTrials, nBlocks])';
    data{iSub}.outcomes = reshape(rawData{iSub}.outcomes, [nTrials, nBlocks])';
    data{iSub}.randHC = reshape(rawData{iSub}.randHC, [nTrials, nBlocks])';
    data{iSub}.randLC = reshape(rawData{iSub}.randLC, [nTrials, nBlocks])';
    data{iSub}.randomReward = reshape(rawData{iSub}.randomReward, [nTrials, nBlocks])';
    data{iSub}.controllability = reshape(rawData{iSub}.controllability, [nTrials, nBlocks])';
end

% ----------------------------------------------------------------------- %
%% Save as one big rawData structure:

fprintf('Save data\n')
outputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
save(outputFile, 'data');

% END OF FILE