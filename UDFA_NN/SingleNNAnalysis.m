function SingleNNAnalysis(loadPath)

    %close all

    % Control parameters
    controlparameters.fs = 48e3;
    controlparameters.nfft = 8192;
    controlparameters.filterOrder = 2;
    controlparameters.nBands = 8;
    controlparameters.fvec = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    controlparameters.fidx = CreateFidx(controlparameters.fvec, controlparameters.nBands);
    epochSize = 20e3;
    
    load([cd filesep loadPath], 'net', 'losses', 'nP');

    disp(['Loss: ', num2str(losses.test(end))])

    controlparameters.numNNInputs = nP.numInputs;
    [inputData, targetData, validationData, geometry] = CreateUDFA_NNTrainingData(epochSize, controlparameters, true, 'ValidationData');
    
    X = dlarray(single(inputData), "CB");

    tfmag = MakeUDFA_NNPrediction(net, X, controlparameters);

    tfmagN = CreateFrequencyNBands(tfmag, controlparameters.fvec, controlparameters.nBands);

    losses.all = CalculateLoss(targetData, tfmagN);

    idx = geometry.zA > 0 & geometry.zA < geometry.wedgeLength;
    percentiles = 10:10:50;

    figure
    PlotVarDependentLoss(losses.all.i, geometry.zA, 1, percentiles, 'zA', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.wedgeLength, 1, percentiles, 'zW', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.bendingAngle, 1, percentiles, 'bA', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.minAngle, 1, percentiles, 'mA', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.rS, 1, percentiles, 'rS', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.rR, 1, percentiles, 'rR', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.zS, 1, percentiles, 'zS', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.zR, 1, percentiles, 'zR', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.thetaS, 1, percentiles, 'tS', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.thetaR, 1, percentiles, 'tR', 1)

    figure
    PlotVarDependentLoss(losses.all.i, geometry.wedgeIndex, 1, percentiles, 'tW', 1)

    figure
    PlotVarDependentLoss(losses.all.i, sum(inputData(1:2,:), 1), 1, percentiles, 'freq 1', 1)

    figure
    PlotVarDependentLoss(losses.all.i, sum(inputData(3:4,:), 1), 1, percentiles, 'freq 2', 1)

    figure
    PlotVarDependentLoss(losses.all.i, sum(inputData(5:6,:), 1), 1, percentiles, 'gains 1', 1)

    figure
    PlotVarDependentLoss(losses.all.i, sum(inputData(7:8,:), 1), 1, percentiles, 'gains 2', 1)
end