classdef LearningModel
    properties
        Q
        epsilon
        beta
        rho
    end

    methods
        function obj = LearningModel(epsilon, rho, beta, numStates, numActions)
            obj.epsilon =epsilon;
            obj.rho = rho;
            obj.beta = beta;
            obj.Q = zeros