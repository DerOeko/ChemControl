function p1 = stableSoftmax(g, ng)
    % Find the maximum of g and ng for numerical stability
    m = max(g, ng);

    g = g - m;
    ng = ng - m;
    
    % Compute the stable softmax components
    exp_g = exp(g);
    exp_ng = exp(ng);
    
    % Calculate the stable softmax output
    p1 = exp_g / (exp_g + exp_ng);
end
