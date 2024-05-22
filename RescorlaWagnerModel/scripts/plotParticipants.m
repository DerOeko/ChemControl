% Clear previous sessions for a clean start
close all;
clear all;
lineStyles = {'-', '--', ':', '-.'}; % Plotting line styles

%% Data Key:
% states: 1 - 8 (5-8 for low control, 1-4 for high control), int8
% controllArray: logical (1=high)
% actions: "None" or "space", logical (1 = Go)
% rewards: 10, 0, -10, int8
% rewardsBool: 1 = reward given/loss avoided in high control trials when correct response is given; in low control trials regardless response, 
% 0 = no reward/loss in high control trials when correct response is given; in low control trials regardless of response
% logical

%% Set up the Import Options and import the data
% Import the data
run("github_config.m") % this includes folder paths with real data, redacted for privacy security
% contains the folderPath
% For this to work, create a github_config.m file with your path to data
% (see README)

fileList = dir(fullfile(folderPath, "*.csv"));
participantData = struct();

for i = 1:numel(fileList)
    fileName = fileList(i).name;
    [~, name, ~] = fileparts(fileName);
    participantId = string(erase(name(1:7), "-"));

    filePath = fullfile(folderPath, fileName);
    opts = detectImportOptions(filePath); % Automatically detect options
    % Specify range and delimiter and variables
    opts.DataLines = [85, Inf];
    opts.Delimiter = ",";
    opts.MissingRule = "omitrow";
    opts.SelectedVariableNames = ["trialType", "miniBlock", "controllability", "randomReward", "feedback", "key_resp_keys"];
    data = readtable(filePath, opts);
    data = renamevars(data, ["trialType", "key_resp_keys"], ["state", "action"]);
    data.action = data.action == "space";
    data.state = int8(data.state);
    data.randomReward = data.randomReward == 1;
    data.controllability = data.controllability == "high";
    data.feedback = int8(data.feedback);
    data.miniBlock = int8(data.miniBlock);
    
    LCdata = data(~data.controllability, :);
    HCdata = data(data.controllability, :);
    participantData.(participantId).LCdata = LCdata;
    participantData.(participantId).HCdata = HCdata;
end

% Clear temporary variables
clear opts d

participants = fieldnames(participantData);
numStates = 4;
proportionGoResponses = zeros(numStates, 10, numel(participants));

% HC Data
for i = 1:numel(participants)
    sub = string(participants(i));
    d = participantData.(sub).HCdata;
    for s = 1:numStates
            stateData = d(d.state == s, :);
            uniqueMiniblocks = unique(stateData.miniBlock);
            numMiniblocks = numel(uniqueMiniblocks);
            for m = 1:numMiniblocks
                miniBlockData = stateData(stateData.miniBlock == uniqueMiniblocks(m, 1), :);
                d.stateOccurrence(d.state == s & d.miniBlock == uniqueMiniblocks(m, 1)) = (1:10)';
            end

            for occ = 1:10
                proportionGoResponses(s, occ, i) = size(d(d.state == s & d.action & d.stateOccurrence == occ, :), 1)/numMiniblocks;
            end
    end
    participantData.(sub).HCdata = d;
end

averageProportionGoResponses = mean(proportionGoResponses, 3);
T = 40;
figure;

hold on; % Allows multiple plots on the same figure
for state = 1:4
    plot(1:T/4, averageProportionGoResponses(state, :)', 'LineWidth', 2); % Plotting mean probabilities for each state
end
xlabel('State Repetitions');
xlim([1.0 T/4])
ylabel('P(Go response | state)');
ylim([0.0, 1.0])
yline(0.5, ":", 'LineWidth', 3, 'Color', '#AEAEAE')
legend('Go to Win', 'Go to Avoid Loss', 'NoGo to Win', 'NoGo to Avoid Loss', 'Location', 'best');
title_str = sprintf("Participant Data: \nMean P(Go|State) Across State Repetitions in\nHigh Control Trials with N = %i", numel(participants));
title(title_str);
grid on

hold off;


