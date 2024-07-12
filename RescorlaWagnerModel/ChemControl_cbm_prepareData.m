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
nBlocks = 9; % how many blocks in one experiment
nTrials = 40; % how many trials in a block

% Downstream settings:

% ----------------------------------------------------------------------- %
%% Retrieve subject IDs:

subFind = dir(fullfile(dirs.behav, "*.csv"));
subFind = {subFind.name};

invalidSubs = [2 10 19 20 38 39 40];
validSubs = setdiff(1:numel(subFind), invalidSubs);
subFind = {subFind{validSubs}};
nSub = numel(subFind);

subList = nan(numel(subFind), 1);
% Retrieve names of subject fiels:
for iSub = 1:nSub % 1 till number of elements in subFind
    subList(iSub) = str2double(string(extractBetween(subFind{iSub}, 'sub-', '_Robo'))); 
    % extracts the actual "sub" number of each subject --> check which ones found
end


% 2, 12, 21, and 22 are below 55% accuracy in high control trials
% 2, 11, 19, 20
% ----------------------------------------------------------------------- %
%% Load rawData:

% Preallocate a cell array to store all the data
rawData = cell(nSub, 1);

% Load raw data
fprintf("Load raw data\n");

for iSub = 1:nSub
    fprintf('Load data for subject %03d\n', iSub);
    fprintf('SubID: %03d\n', subList(iSub));
    % Define file path
    filePath = fullfile(folderPath, subFind{iSub});
    
    % Define options for reading calibration data (rows 44 to 83)
    optsCal = detectImportOptions(filePath);
    optsCal.DataLines = [44, 83]; % Lines 44 to 83
    optsCal.Delimiter = ","; % CSV delimiter
    optsCal.MissingRule = "fill"; % Fill missing data
    optsCal.SelectedVariableNames = ["trialType", "miniBlock", "controllability", "randHC", "randLC", "randomReward", "feedback", "key_respCalibrate_keys"];
    
    % Define options for reading main data (rows 85 to Inf)
    optsMain = detectImportOptions(filePath);
    optsMain.DataLines = [85, Inf]; % Start reading from line 85
    optsMain.Delimiter = ","; % CSV delimiter
    optsMain.MissingRule = "omitrow"; % Omit rows with missing data
    optsMain.SelectedVariableNames = ["trialType", "miniBlock", "controllability", "randHC", "randLC", "randomReward", "feedback", "key_resp_keys"];
    
    % Try to load the calibration data
    try
        calibrationData = readtable(filePath, optsCal);
    catch ME
        fprintf('Error loading calibration data for subject %03d: %s\n', iSub, ME.message);
        calibrationData = table(); % Create an empty table if an error occurs
    end
    
    % Fill missing data in calibrationData
    calibrationData.randLC(:) = 2;

    % Rename variable names
    calibrationData = renamevars(calibrationData, ["trialType", "key_respCalibrate_keys", "feedback"], ["stimuli", "actions", "outcomes"]);

    % Try to load the main data
    try
        mainData = readtable(filePath, optsMain);
    catch ME
        fprintf('Error loading main data for subject %03d: %s\n', iSub, ME.message);
        mainData = table(); % Create an empty table if an error occurs
    end

    mainData = renamevars(mainData, ["trialType", "key_resp_keys", "feedback"], ["stimuli", "actions", "outcomes"]);

    % Combine calibration data and main data
    data = [calibrationData; mainData];
    
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
    
    % Convert outcomes to -1, 0, 1
    data.outcomes = double((data.outcomes == 10) - (data.outcomes == -10));
    
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