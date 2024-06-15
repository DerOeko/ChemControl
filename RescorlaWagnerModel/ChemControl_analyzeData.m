% ChemControl_analyzeData.m
% 
% Execute this script to run analyses for the preparedData.
% ChemControl_cbm_prepareData.m has to be run before.
%
% INPUTS:
% none.
%
% OUTPUTS:
% No outputs, just plots
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/

% clear all; close all; clc

%% Directories

dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

%% Input file

inputFile = fullfile(dirs.target, 'ChemControl_cbm_inputData.mat');
data = load(inputFile).data;
fprintf('Loaded inputFile with %d subjects\n', size(data, 2));

%% Settings
nSub = size(data,2);
nBlocks = size(data{1}.actions, 1);
nTrials = size(data{1}.actions, 2);

%% PLOTS
        


