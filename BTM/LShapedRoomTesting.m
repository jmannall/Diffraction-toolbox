close all    
clear all

    x = [8 14];
    y = [6 10];
    height = 3;

    fs = 48e3;
    nfft = 4096;
    difforder = 1;
    createPlot = false;
    % receiver = [sum(x) sum(y) height] / 2;

    numX = 200;
    numY = 200;
    numReceivers = numX * numY;
    epsilon = 1e-5;
    distX = linspace(epsilon, x(2) - epsilon, numX)';
    distY = linspace(epsilon, y(2) - epsilon, numY)';
    
    count = 0;
    receiver = zeros(numReceivers, 3);
    for i = 1:numX
        for j = 1:numY
            count = count + 1;
            receiver(count,:) = [distX(i) distY(j) height / 2];
        end
    end

    outRoom = receiver(:,1) >= x(1) & receiver(:,2) <= y(1);
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', difforder);


    %[tfmag, fvec, tfcomplex] = SingleLShapedRoom(x, y, height, receiver, controlparameters, createPlot);

    parfor i = 1:numReceivers
        [tfmag(:,i), fvec(:,i), tfcomplex(:,i)] = SingleLShapedRoom(x, y, height, receiver(i,:), controlparameters, outRoom(i), createPlot);
    end

    idx = find(fvec(:,1) > 300, 1, "first");
    tfcomplexAll = [tfcomplex.complete];
    tfcomplexDir = [tfcomplex.direct];
    tfcomplexSpec = [tfcomplex.geom];
    tfcomplexDiff = [tfcomplex.diff1];

    tfcomplexAll = reshape(tfcomplexAll, nfft / 2, [], numX);
    tfcomplexDir = reshape(tfcomplexDir, nfft / 2, [], numX);
    tfcomplexSpec = reshape(tfcomplexSpec, nfft / 2, [], numX);
    tfcomplexDiff = reshape(tfcomplexDiff, nfft / 2, [], numX);

    PlotSpectrogram(squeeze(tfcomplexAll(idx,:,:)), distX, distY, [-70 0], 'All 300Hz', false, false)
    PlotSpectrogram(squeeze(tfcomplexDir(idx,:,:)), distX, distY, [-70 0], 'Dir 300Hz', false, false)
    PlotSpectrogram(squeeze(tfcomplexSpec(idx,:,:)), distX, distY, [-70 0], 'Spec 300Hz', false, false)
    PlotSpectrogram(squeeze(tfcomplexDiff(idx,:,:)), distX, distY, [-70 0], 'Diff 300Hz', false, false)