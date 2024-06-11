dirs.root           = '/project/3017083.01/behavioral_study/scripts/matlab_scripts/RescorlaWagnerModel';
dirs.models         = fullfile(dirs.root, 'models');
nMod = length(dir(fullfile(dirs.models, '*.m')));

% fig1 = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

% for iMod = 1:nMod
%     subplot(3, 3, iMod);
%     subj = sim_subj();
%     parameters = [0.05 1.5 0.5 0.3 0.2 1.5 1.5];
%     output = eval(sprintf("ChemControl_mod%d_modSim(parameters, subj)", iMod));
%     HCcell = output.HCcell;
%     occurrenceMeans = retrieveOccurrenceMeans(HCcell);
%     plotLearningCurves(occurrenceMeans, sprintf("ChemControl_mod%d_modSim", iMod), true, fig1);
% end
T = 40;
R = 50;

fprintf('Initialize directories\n');
dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
fdata = load(inputFile);
data = fdata.data;
nSub = length(data);

HCmeans = zeros(T/4, 4, nSub);
parameters = [0.07, 17, 0.45, 0.75];

for iSub = 1:nSub
    subj = data{nSub};
    out = ChemControl_mod3_modSim(parameters, subj);
    HCcell = out.HCcell;

    HCoccurrences = zeros(T/4, 4);
    for s = 1:4
        for occurrence = 1:T/4
            HCprobs = zeros(B/2, 1);
            for b = 1:B/2
                HCprobs(b) = HCcell{b, s}(occurrence);
            end
            HCoccurrences(occurrence, s) = mean(HCprobs);
        end
    end
    HCmeans(:, :, iSub) = HCoccurrences;
end

HCmeans = mean(HCmeans, 3);

plotLearningCurves(HCmeans, "FixedPavlov", true);
