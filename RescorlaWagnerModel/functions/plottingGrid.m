function figHandle = plottingGrid(epsilonFlat, betaFlat, rhoFlat, accuracyFlat)

    if nargin < 4
        disp('Missing arguments. Running runGridSearch to generate data...');
        [accuracyMatrix, epsilonFlat, betaFlat, rhoFlat, accuracyFlat] = runGridSearch();
    end
    figure;
    figHandle = figure; % Capture the handle for the created figure

    scatter3(epsilonFlat, betaFlat, rhoFlat, 36, accuracyFlat, "filled");
    title('Grid Search Results');
    xlabel('Epsilon');
    ylabel('Beta');
    zlabel('Rho');
    cb = colorbar; 
    ylabel(cb, 'Accuracy');
    colormap('jet');
    view(3);
end
