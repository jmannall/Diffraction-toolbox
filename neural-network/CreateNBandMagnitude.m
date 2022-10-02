function tfmagNBand = CreateNBandMagnitude(tfmag, fidx)
    numBands = max(fidx);
    numObservations = size(tfmag, 1);
    tfmagNBand = dlarray(zeros(numObservations, numBands));
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(:,i) = sum(tfmag(:,fidx == i), 2) / num;
    end
end