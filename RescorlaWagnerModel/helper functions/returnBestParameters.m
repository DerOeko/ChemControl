function parameters = returnBestParameters(group_means, transforms)
%RETURNBESTPARAMETERS Returns the best parameters of a model by applying
%the transforms specified on the group means.
parameters = zeros(length(group_means), 1);

for i = 1:length(group_means)
    param_value = group_means(i);
    transform_fn = str2func(transforms{i});
    transformed_value = transform_fn(param_value);
    parameters(i) = transformed_value;
end

end

