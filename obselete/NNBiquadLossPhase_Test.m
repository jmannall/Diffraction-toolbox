function [loss, lossAll] = NNBiquadLossPhase_Test(net, inputData, targetData, phaseData, numBiquads, numFreq, numObservations, fidx)
    
    targetMag = targetData;
    targetPhase = phaseData;
    numObservations = round(numObservations);
    output = predict(net, inputData);
    
    [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
    [tfmag, fvec, tfcomplex] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

    phase = angle(tfcomplex);
    c = 343;
    d = 2;
    delay = d / c;
    phaseShift = -2 * pi * delay * fvec;
    phase = phase + phaseShift;
    phaseLoss = mean((phase - targetPhase).^2, "all");

    numBands = max(fidx);
    tfmagNBand = dlarray(zeros(numObservations, numBands));
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(:,i) = sum(tfmag(:,fidx == i), 2) / num;
    end

    %loss = (tfmagNBand' - 128 * targetData).^2 / (numBands * numObservations);
    %loss = sum(loss, "all");
    lossAll = sqrt(mean((tfmagNBand' - 128 * targetMag).^2));
    loss = sqrt(mean((tfmagNBand' - 128 * targetMag).^2, "all"));
end