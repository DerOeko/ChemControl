function [subj] = sim_subj()

% sim_subj

% 
% Initialize data for one "standard" subject for model simulations.
%
% INPUTS:
% none.
%
% OUTPUTS:
% subj          = structure with the following fields:
% .stimuli      = vector of integers, stimuli seen (1-16).
% .outcomes 	= vector of integers, 
%
% CHEMCONTROL STUDY, DONDERS INSTITUTE, NIJMEGEN.
% S. Nellessen, 2024.
% Inspired by code from Algermissen, J. et al. 2024.

% we are here:
% cd dirs.root/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% Settings:
B = 8; % number of blocks
T = 40; % number of trials in a block
cP = 0.8; % control probability
rP = 0.85; % reward probability
numControl = round(cP * T); % number of trials with control
numRewarded = round(rP * T); % number of rewarded trials

% Cues and how often they repeat
cond = [1, 2, 3, 4]; 
rep = T / 4;

% Initialize randLC, randHC, controllabilities, stimuli, outcomes, actions,
% randomrewards, 
stimuli = zeros(B, T);
randHC = zeros(B, T);
randLC = zeros(B, T);
controllability = zeros(B, T);
randomReward = zeros(B, T);
cArr = repelem([1 0], B/2); 
cArr = cArr(randperm(B));
for b = 1:B
    sVec = repelem(cond, rep);
    sVec = sVec(randperm(length(stimuli)));
    stimuli(b, :) = sVec;
    
    rVec = [ones(1, numRewarded) zeros(1, T-numRewarded, 1)];
    rVec = rVec(randperm(T));
    randomReward(b, :) = rVec;
    if cArr(b)
        controllability(b, :) = 1;
        cVec = [ones(1, numControl), zeros(1, T-numControl, 1)];
        cVec = cVec(randperm(T));

        randHC(b, :) = cVec;
        randLC(b, :) = 2;
    else
        controllability(b, :) = 0;
        cVec = [ones(1, numControl), zeros(1, T-numControl, 1)];
        cVec = cVec(randperm(T));

        randLC(b, :) = cVec;
        randHC(b, :) = 2;
    end
end

subj.stimuli = stimuli;
subj.randHC = randHC;
subj.randLC = randLC;
subj.controllability = controllability;
subj.randomReward = randomReward;

end




    
    


