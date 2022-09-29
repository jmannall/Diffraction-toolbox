function [loss, lossAll] = NNBiquadLoss_Test(net, inputData, targetData, numBiquads, numFreq, numObservations, fidx)

    numObservations = round(numObservations);
    output = predict(net, inputData);
    
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
    lossAll = sqrt(mean((tfmagNBand' - 128 * targetData).^2));
    loss = sqrt(mean((tfmagNBand' - 128 * targetData).^2, "all"));
end