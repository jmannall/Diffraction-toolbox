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
        tfmagNBand(i,:) = mag2db(abs(sum(db2mag(tfmag(fidx == i,:)), 1) / num));
        %tfmagNBand(i,:) = sum(tfmag(fidx == i,:), 1) / num;     %
        %incorrect need to average as magnitude not dB
    end
end