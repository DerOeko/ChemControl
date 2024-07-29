function [hc_d,lc_d,yc_d] = extractControlTypeData(data)
nSub = numel(data);
%% Extracting subsets of data
hc_d = {};
lc_d = {};
yc_d = {};

for iSub = 1:nSub
    d = data{iSub};  % Get the struct for the current subject
    
    % Identify high control (hc) blocks
    hc_idx = find(d.controllability(:, 1));
    
    % Identify low control (lc) blocks where controllability is 0 and isYoked is 0
    lc_idx = find(~d.controllability(:, 1) & ~d.isYoked(:, 1));
    
    % Identify yoked control (yc) blocks where isYoked is 1
    yc_idx = find(d.isYoked(:, 1));
    
    % Extract the data for the high control blocks and store it
    hc_d{iSub} = structfun(@(x) x(hc_idx, :), d, 'UniformOutput', false);
    lc_d{iSub} = structfun(@(x) x(lc_idx, :), d, 'UniformOutput', false);
    yc_d{iSub} = structfun(@(x) x(yc_idx, :), d, 'UniformOutput', false);
end

end

