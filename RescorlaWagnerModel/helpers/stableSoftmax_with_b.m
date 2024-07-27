function p1 = stableSoftmax_with_b(g, ng, beta)
    % Find the maximum of g and ng for numerical stability
    m = max(g, ng);
    g = g - m;
    ng = ng - m;
    
    % Compute the stable softmax components
    exp_g = exp(g*beta);
    exp_ng = exp(ng*beta);
    
    % Calculate the stable softmax output
    p1 = (exp_g) / (exp_g + exp_ng);
end
