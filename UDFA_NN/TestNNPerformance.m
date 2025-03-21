function TestNNPerformance(numData)

    %close all

    % Control parameters
    controlparameters.fs = 48e3;
    controlparameters.nfft = 8192;
    controlparameters.filterOrder = 2;
    controlparameters.nBands = 8;
    controlparameters.fvec = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    controlparameters.fidx = CreateFidx(controlparameters.fvec, controlparameters.nBands);

    % Load paths
    rootDir = 'NNSaves';
    loadDir = 'SingleUDFA_NN';
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep loadDir])
    
    loadPath = [rootDir filesep loadDir];
    files = dir([cd filesep loadPath]);
    files = files(~[files.isdir]);

    numFiles = length(files);

    colours = colororder;
    colours = colours(mod(0:numData - 1, length(colours)) + 1,:);
    for i = 1:numFiles
        if (numFiles == 1)
            file = files;
        else
            file = files(i);
        end
        load([cd filesep loadPath filesep file.name], 'net', 'losses', 'nP');

        disp(file.name)
        disp(['Loss: ', num2str(losses.test(end))])

        controlparameters.numNNInputs = nP.numInputs;
        [inputData, ~, validationData] = CreateSingleUDFA_NNTrainingData(numData, controlparameters, true, 'ValidationData');
        X = dlarray(single(inputData), "CB");

        tfmag = MakeUDFA_NNPrediction(net, X, controlparameters);

        PlotNNTrainingLossess(losses.iteration, losses.epoch, losses.test, file.name)

        figure
        colororder(colours);
        semilogx(controlparameters.fvec, validationData)
        hold on
        grid on
        semilogx(controlparameters.fvec, tfmag, '--')
        title(replace(file.name, '_', ' '))
        ylim([-70 10])
        xlim([20 20e3])
    end
end