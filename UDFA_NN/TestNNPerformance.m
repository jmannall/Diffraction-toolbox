function TestNNPerformance(loadDir)

    close all

    % Control parameters
    controlparameters.fs = 48e3;
    controlparameters.nfft = 8192;
    controlparameters.filterOrder = 2;
    controlparameters.nBands = 8;
    controlparameters.fvec = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    controlparameters.fidx = CreateFidx(controlparameters.fvec, controlparameters.nBands);

    % Load paths
    rootDir = 'NNSaves';
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep loadDir])
    
    loadPath = [rootDir filesep loadDir];
    files = dir([cd filesep loadPath]);
    files = files(~[files.isdir]);

    numFiles = length(files);
    numData = 5;

    colours = colororder;
    colours = colours(1:numData,:);
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
        [inputData, ~, validationData] = CreateUDFATrainingData(5, controlparameters, 'UDFA_NN', true, 'ValidationData');
        X = dlarray(single(inputData), "CB");

        tfmag = MakeUDFA_NNPrediction(net, X, controlparameters);

        PlotNNTrainingLossess(losses.iteration, losses.epoch, file.name)

        figure
        colororder(colours);
        semilogx(controlparameters.fvec, validationData)
        hold on
        grid on
        semilogx(controlparameters.fvec, tfmag, '--')
        title(file.name)
    end
end