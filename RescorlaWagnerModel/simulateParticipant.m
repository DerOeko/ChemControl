% Clear previous sessions for a clean start
close all;
clear all;

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 228, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = [85, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "trialType", "Var4", "Var5", "controllability", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61", "Var62", "Var63", "Var64", "Var65", "Var66", "Var67", "Var68", "Var69", "Var70", "Var71", "Var72", "Var73", "Var74", "Var75", "Var76", "Var77", "Var78", "Var79", "Var80", "Var81", "Var82", "Var83", "randomReward", "feedback1", "Var86", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "Var93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102", "Var103", "Var104", "Var105", "Var106", "Var107", "Var108", "Var109", "Var110", "Var111", "Var112", "Var113", "Var114", "Var115", "Var116", "Var117", "Var118", "Var119", "Var120", "Var121", "Var122", "Var123", "Var124", "Var125", "Var126", "Var127", "Var128", "Var129", "Var130", "Var131", "Var132", "Var133", "Var134", "Var135", "Var136", "Var137", "Var138", "Var139", "Var140", "Var141", "Var142", "Var143", "Var144", "Var145", "Var146", "Var147", "Var148", "Var149", "Var150", "Var151", "Var152", "Var153", "Var154", "Var155", "Var156", "Var157", "Var158", "Var159", "Var160", "Var161", "Var162", "Var163", "Var164", "Var165", "Var166", "Var167", "Var168", "Var169", "Var170", "Var171", "Var172", "Var173", "Var174", "Var175", "Var176", "Var177", "Var178", "Var179", "Var180", "Var181", "Var182", "Var183", "Var184", "Var185", "Var186", "Var187", "Var188", "key_respkeys", "Var190", "Var191", "Var192", "Var193", "Var194", "Var195", "Var196", "Var197", "Var198", "Var199", "Var200", "Var201", "Var202", "Var203", "Var204", "Var205", "Var206", "Var207", "Var208", "Var209", "Var210", "Var211", "Var212", "Var213", "Var214", "Var215", "Var216", "Var217", "Var218", "Var219", "Var220", "Var221", "Var222", "Var223", "Var224", "Var225", "Var226", "Var227", "Var228"];
opts.SelectedVariableNames = ["trialType", "controllability", "randomReward", "feedback1", "key_respkeys"];
opts.VariableTypes = ["string", "string", "int8", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "logical", "int8", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "categorical", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var5", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61", "Var62", "Var63", "Var64", "Var65", "Var66", "Var67", "Var68", "Var69", "Var70", "Var71", "Var72", "Var73", "Var74", "Var75", "Var76", "Var77", "Var78", "Var79", "Var80", "Var81", "Var82", "Var83", "Var86", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "Var93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102", "Var103", "Var104", "Var105", "Var106", "Var107", "Var108", "Var109", "Var110", "Var111", "Var112", "Var113", "Var114", "Var115", "Var116", "Var117", "Var118", "Var119", "Var120", "Var121", "Var122", "Var123", "Var124", "Var125", "Var126", "Var127", "Var128", "Var129", "Var130", "Var131", "Var132", "Var133", "Var134", "Var135", "Var136", "Var137", "Var138", "Var139", "Var140", "Var141", "Var142", "Var143", "Var144", "Var145", "Var146", "Var147", "Var148", "Var149", "Var150", "Var151", "Var152", "Var153", "Var154", "Var155", "Var156", "Var157", "Var158", "Var159", "Var160", "Var161", "Var162", "Var163", "Var164", "Var165", "Var166", "Var167", "Var168", "Var169", "Var170", "Var171", "Var172", "Var173", "Var174", "Var175", "Var176", "Var177", "Var178", "Var179", "Var180", "Var181", "Var182", "Var183", "Var184", "Var185", "Var186", "Var187", "Var188", "Var190", "Var191", "Var192", "Var193", "Var194", "Var195", "Var196", "Var197", "Var198", "Var199", "Var200", "Var201", "Var202", "Var203", "Var204", "Var205", "Var206", "Var207", "Var208", "Var209", "Var210", "Var211", "Var212", "Var213", "Var214", "Var215", "Var216", "Var217", "Var218", "Var219", "Var220", "Var221", "Var222", "Var223", "Var224", "Var225", "Var226", "Var227", "Var228"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var5", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61", "Var62", "Var63", "Var64", "Var65", "Var66", "Var67", "Var68", "Var69", "Var70", "Var71", "Var72", "Var73", "Var74", "Var75", "Var76", "Var77", "Var78", "Var79", "Var80", "Var81", "Var82", "Var83", "Var86", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "Var93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102", "Var103", "Var104", "Var105", "Var106", "Var107", "Var108", "Var109", "Var110", "Var111", "Var112", "Var113", "Var114", "Var115", "Var116", "Var117", "Var118", "Var119", "Var120", "Var121", "Var122", "Var123", "Var124", "Var125", "Var126", "Var127", "Var128", "Var129", "Var130", "Var131", "Var132", "Var133", "Var134", "Var135", "Var136", "Var137", "Var138", "Var139", "Var140", "Var141", "Var142", "Var143", "Var144", "Var145", "Var146", "Var147", "Var148", "Var149", "Var150", "Var151", "Var152", "Var153", "Var154", "Var155", "Var156", "Var157", "Var158", "Var159", "Var160", "Var161", "Var162", "Var163", "Var164", "Var165", "Var166", "Var167", "Var168", "Var169", "Var170", "Var171", "Var172", "Var173", "Var174", "Var175", "Var176", "Var177", "Var178", "Var179", "Var180", "Var181", "Var182", "Var183", "Var184", "Var185", "Var186", "Var187", "Var188", "key_respkeys", "Var190", "Var191", "Var192", "Var193", "Var194", "Var195", "Var196", "Var197", "Var198", "Var199", "Var200", "Var201", "Var202", "Var203", "Var204", "Var205", "Var206", "Var207", "Var208", "Var209", "Var210", "Var211", "Var212", "Var213", "Var214", "Var215", "Var216", "Var217", "Var218", "Var219", "Var220", "Var221", "Var222", "Var223", "Var224", "Var225", "Var226", "Var227", "Var228"], "EmptyFieldRule", "auto");

% Import the data
d = readtable("/project/3017083.01/behavioral study/data/raw/sub-X11_RoboMadness_Behave_2024-03-01_15h12.44.602.csv", opts);

%% Convert to output type
% the stimuli in each action-state pair evaluation the model learns
% 2 valences * 2 actions * 2 controllability conditions = 8
states = d.trialType;
states = rmmissing(states);

% 1 = reward given/loss avoided in high control trials when correct response is given; in low control trials regardless response, 
% 0 = no reward/loss in high control trials when correct response is given; in low control trials regardless of response
rewardsBool = d.randomReward;
rewardsBool = rmmissing(rewardsBool);

% the actual reward, to calculate the prediction error
rewards = d.feedback1;
rewards = rmmissing(rewards);

% space or None, depending on the response of the participant
actions = d.key_respkeys;
actions = rmmissing(actions);

% High or Low controllability
controllArray = d.controllability;
controllArray = rmmissing(controllArray);

% Clear temporary variables
clear opts d

%% Data Key:
% states: 1 - 8 (5-8 for low control, 1-4 for high control), int8
% controllArray: "high" or "low", string
% actions: "None" or "space", categorical
% rewards: 10, 0, -10, int8
% rewardsBool: 1 = reward given/loss avoided in high control trials when correct response is given; in low control trials regardless response, 
% 0 = no reward/loss in high control trials when correct response is given; in low control trials regardless of response
% logical

%% split high and low control trials
disp(controllArray)
HCfilterArray = strcmp(controllArray, "high");

HCstates = states(HCfilterArray);
HCactions = actions(HCfilterArray);
HCrewards = rewards(HCfilterArray);
HCrewardsBool = rewardsBool(HCfilterArray);

LCstates = states(~HCfilterArray);
LCactions = actions(~HCfilterArray);
LCrewards = rewards(~HCfilterArray);
LCrewardsBool = rewardsBool(~HCfilterArray);
