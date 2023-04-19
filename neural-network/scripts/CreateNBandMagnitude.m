function tfmagNBand = CreateNBandMagnitude(tfmag, fidx)
    numBands = max(fidx);
    numObservations = size(tfmag, 2);
%     if isdlarray(tfmag)
%         tfmagNBand = dlarray(zeros(numBands, numObservations));
%     else
%         tfmagNBand = zeros(numBands, numObservations);
%     end
    tfmagNBand = zeros(numBands, numObservations);
    for i = 1:numBands
        num = sum(fidx == i);
        tfmagNBand(i,:) = 20*log(abs(sum(10 .^ (tfmag(fidx == i,:) / 20), 1) / num)) / log(10);
        %tfmagNBand(i,:) = sum(tfmag(fidx == i,:), 1) / num;     %
        %incorrect need to average as magnitude not dB
    end
end