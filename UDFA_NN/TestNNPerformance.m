function TestNNPerformance(numData)
%close all
set(gca,'ColorOrder','factory')
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
    loadDir = ['UDFA_NN' filesep 'Run0'];
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep loadDir])
    
    loadPath = [rootDir filesep loadDir];
    files = dir([cd filesep loadPath]);
    files = files(~[files.isdir]);

    numFiles = length(files);

    colours = colororder;
    colours = colours(mod(0:numData - 1, length(colours)) + 1,:);
    for i = 5:5
        if (numFiles == 1)
            file = files;
        else
            file = files(i);
        end
        load([cd filesep loadPath filesep file.name], 'net', 'losses', 'nP');

        disp(file.name)
        disp(['Loss: ', num2str(losses.test(end))])

        controlparameters.numNNInputs = nP.numInputs;
        [inputData, ~, validationData] = CreateUDFA_NNTrainingData(numData, controlparameters, true, 'ValidationData');
        X = dlarray(single(inputData), "CB");

        tfmag = MakeUDFA_NNPrediction(net, X, controlparameters);

        %PlotNNTrainingLossess(losses.iteration, losses.epoch, losses.test, file.name)

        for j = 1:numData
            [b, a] = fracOrderBlendLPapprox3(10 ^ inputData(j), 0.5, 1.44, 0.2, 2, controlparameters.fs);
            udfaIIR2(:,j) = CalculateFilterResponse(b', a', controlparameters.nfft, controlparameters.fs);
            [b, a] = fracOrderBlendLPapprox3(10 ^ inputData(j), 0.5, 1.44, 0.2, 4, controlparameters.fs);
            udfaIIR4(:,j) = CalculateFilterResponse(b', a', controlparameters.nfft, controlparameters.fs);
        end

        figure
        colororder(colours);
        semilogx(controlparameters.fvec, validationData, '-k', 'LineWidth', 3)
        hold on
        grid on
        semilogx(controlparameters.fvec, tfmag, '-r')
        semilogx(controlparameters.fvec, udfaIIR2, '--b')
        semilogx(controlparameters.fvec, udfaIIR4, '-.g')
        title(replace(file.name, '_', ' '))
        ylim([-40 10])
        xlim([20 20e3])

        target = CreateNBandMagnitude(validationData, controlparameters.fidx);
        tfmagNBand = CreateNBandMagnitude(tfmag, controlparameters.fidx);
        lossNN(i) = mean((tfmagNBand - target).^2, 'all');
        tfmagNBand = CreateNBandMagnitude(udfaIIR2, controlparameters.fidx);
        lossUDFA2(i) = mean((tfmagNBand - target).^2, 'all');
        tfmagNBand = CreateNBandMagnitude(udfaIIR4, controlparameters.fidx);
        lossUDFA4(i) = mean((tfmagNBand - target).^2, 'all');
    end
end