function TestNNPerformance(loadPath, numData)

    %close all

    % Control parameters
    controlparameters.fs = 48e3;
    controlparameters.nfft = 8192;
    controlparameters.filterOrder = 2;
    controlparameters.nBands = 8;
    controlparameters.fvec = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    controlparameters.fidx = CreateFidx(controlparameters.fvec, controlparameters.nBands);
    
    load([cd filesep loadPath], 'net', 'losses', 'nP');

    numEpochs = length(losses.epoch);
    numIterations = length(losses.iteration) / numEpochs;
    losses.test = losses.test(losses.test ~= 1000);
    losses.epoch = losses.epoch(losses.epoch ~= 1000);
    losses.iteration = losses.iteration(1:length(losses.epoch) * numIterations / numEpochs);

    disp(['Loss: ', num2str(losses.test(end))])

    controlparameters.numNNInputs = nP.numInputs;
    [inputData, ~, validationData, geometry] = CreateUDFA_NNTrainingData(numData, controlparameters, true, 'ValidationData');
    numData = length(geometry.wedgeIndex);

    X = dlarray(single(inputData), "CB");

    tfmag = MakeUDFA_NNPrediction(net, X, controlparameters);

    PlotNNTrainingLossess(losses.iteration, losses.epoch, losses.test, loadPath)

    colours = colororder;
    colours = colours(mod(0:numData - 1, length(colours)) + 1,:);
    
    figure
    colororder(colours);
    semilogx(controlparameters.fvec, validationData)
    hold on
    grid on
    semilogx(controlparameters.fvec, tfmag, '--')
    legendEntries = [{'Target'}, repmat({''}, 1, numData - 1), {'NN-UDFA'}];
    legend(legendEntries)
    title(replace(loadPath, '_', ' '))
    ylim([-70 10])
    xlim([20 20e3])
end