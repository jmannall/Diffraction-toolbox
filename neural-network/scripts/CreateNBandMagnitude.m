function tfmagNBand = CreateNBandMagnitude(tfmag, fidx)
    numBands = max(fidx);
    numObservations = size(tfmag, 2);
    if isdlarray(tfmag)
        tfmagNBand = dlarray(zeros(numBands, numObservations));
    else
        tfmagNBand = zeros(numBands, numObservations);
    end
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(i,:) = sum(tfmag(fidx == i,:), 1) / num;
    end
end