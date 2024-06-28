function [subj] = sim_subj(B, T)

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
cP = 0.8; % control probability
rP = 0.8; % reward probability
numControl = round(cP * T); % number of trials with control
numRewarded = round(rP * T); % number of rewarded trials
controllabilitySchedules = [
    1, 2, 2, 1, 2, 1, 1, 2;
    2, 1, 2, 2, 1, 2, 1, 1;
    1, 2, 1, 2, 2, 1, 2, 1;
    1, 1, 2, 1, 2, 2, 1, 2;
    2, 1, 1, 2, 1, 2, 2, 1;
    1, 2, 1, 1, 2, 1, 2, 2;
    2, 1, 2, 1, 1, 2, 1, 2;
    2, 2, 1, 2, 1, 1, 2, 1;
    1, 1, 2, 2, 1, 2, 1, 2;
    2, 2, 1, 1, 2, 1, 2, 1;
    1, 2, 1, 2, 1, 2, 1, 2;
];

% Cues and how often they repeat
cond = [1, 2, 3, 4]; 
rep = T / 4;

% Initialize randLC, randHC, controllabilities, stimuli, outcomes, actions,
% randomrewards, 
stimuli = zeros(B, T);
randHC = zeros(B, T);
randLC = zeros(B, T);
randomReward = zeros(B, T);
controllability = zeros(B, T);

sVec = repelem(cond, rep);
sVec = sVec(randperm(size(stimuli, 2)));
cali_stimuli(:) = sVec;
rVec = [ones(1, numRewarded) zeros(1, T-numRewarded, 1)];
rVec = rVec(randperm(T));
cali_randReward(:) = rVec;
cVec = [ones(1, numControl), zeros(1, T-numControl, 1)];
cVec = cVec(randperm(T));
cali_randLC = zeros(1, T);
cali_randHC(:) = cVec;
cali_randLC(:) = 2;

% Select a random controllability schedule
selected_schedule_idx = randi(size(controllabilitySchedules, 1));
cArr = controllabilitySchedules(selected_schedule_idx, :);
two_indices = find(cArr == 2);
selected_indices = randsample(two_indices, 2);
cArr(selected_indices) = 0;

for b = 1:B
    sVec = repelem(cond, rep);
    sVec = sVec(randperm(size(stimuli, 2)));
    stimuli(b, :) = sVec;
    rVec = [ones(1, numRewarded) zeros(1, T-numRewarded, 1)];
    rVec = rVec(randperm(T));
    randomReward(b, :) = rVec;

    index = mod(b-1, length(cArr)) + 1;
    switch cArr(index)
        case 1
            controllability(b, :) = 1;
            cVec = [ones(1, numControl), zeros(1, T-numControl, 1)];
            cVec = cVec(randperm(T));
    
            randHC(b, :) = cVec;
            randLC(b, :) = 2;
        case 0
            controllability(b, :) = 0;
            cVec = [ones(1, numControl), zeros(1, T-numControl, 1)];
            cVec = cVec(randperm(T));
    
            randLC(b, :) = cVec;
            randHC(b, :) = 2;
        case 2
            controllability(b, :) = 2;
            cVec = [ones(1, numControl), zeros(1, T-numControl, 1)];
            cVec = cVec(randperm(T));
    
            randLC(b, :) = cVec;
            randHC(b, :) = 2;
    end
end

subj.stimuli = stimuli;
subj.randHC = randHC;
subj.randLC = randLC;
subj.randomReward = randomReward;
subj.controllability = controllability;
subj.cali_stimuli = cali_stimuli;
subj.cali_randReward = cali_randReward;
subj.cali_randHC = cali_randHC;
subj.cali_randLC = cali_randLC;
subj.selected_schedule_idx = selected_schedule_idx;
subj.selected_schedule = cArr;

end




    
    


