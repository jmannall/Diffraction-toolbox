function [loss, state, gradients] = NNBiquadLossPhase(net, inputData, targetData, phaseData, numBiquads, numFreq, numObservations, fidx)

    targetMag = targetData;
    targetPhase = phaseData;
    [output, state] = forward(net, inputData);
    
    [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
    [tfmag, fvec, tfcomplex] = CreateBiquad(zR, zI, pR, pI, k, numFreq);
    
    phase = angle(tfcomplex);
    magTemp = abs(tfcomplex);
    c = 343;
    d = 2;
    delay = d / c;
    phaseShift = 2 * pi * delay * fvec;
    phase = phase + phaseShift;
    tfcomplexTemp = magTemp .* exp(phase * 1i);
    phase = angle(tfcomplexTemp);
    phaseLoss = mean((phase(1,:) - targetPhase(1,:)).^2, "all");

    numBands = max(fidx);
    tfmagNBand = dlarray(zeros(numObservations, numBands));
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(:,i) = sum(tfmag(:,fidx == i), 2) / num;
    end

    %loss = (tfmagNBand' - 128 * targetData).^2 / (numBands * numObservations);
    %loss = sum(loss, "all");
    magLoss = mean((tfmagNBand' - 128 * targetMag).^2, "all");
    loss = magLoss + 10 * phaseLoss;
    gradients = dlgradient(loss, net.Learnables);
end