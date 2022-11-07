function tfmagNBand = CreateNBandMagnitude(tfmag, fidx)
    numBands = max(fidx);
    numObservations = size(tfmag, 2);
    tfmagNBand = dlarray(zeros(numBands, numObservations));
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(i,:) = sum(tfmag(fidx == i,:), 1) / num;
    end
end