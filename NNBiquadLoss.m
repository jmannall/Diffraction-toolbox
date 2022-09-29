function [loss, state, gradients] = NNBiquadLoss(net, inputData, targetData, numBiquads, numFreq, numObservations, fidx)

    [output, state] = forward(net, inputData);
    
    [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
    [tfmag, ~] = CreateBiquad(zR, zI, pR, pI, k, numFreq);
    
    numBands = max(fidx);
    tfmagNBand = dlarray(zeros(numObservations, numBands));
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(:,i) = sum(tfmag(:,fidx == i), 2) / num;
    end

    %loss = (tfmagNBand' - 128 * targetData).^2 / (numBands * numObservations);
    %loss = sum(loss, "all");
    loss = mean((tfmagNBand' - 128 * targetData).^2, "all");
    gradients = dlgradient(loss, net.Learnables);
end