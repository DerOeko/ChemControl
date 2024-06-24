function p1 = stableSoftmax(g, ng)
    % Find the maximum of g and ng for numerical stability
    m = max(g, ng);
    
    % Compute the stable softmax components
    exp_g = exp(g - m);
    exp_ng = exp(ng - m);
    
    % Calculate the stable softmax output
    p1 = exp_g / (exp_g + exp_ng);
end
