%% IIR filter loss function

function loss = IIRFilterLoss(output, target, numIIRFilters, numFreq, fs, fidx)
    
    [z, p, k] = CreateIIRFromNNOutput(output, numIIRFilters);

    [tfmag, ~] = CreateIIRFilter(z, p, k, numFreq, fs);
    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));
    loss = sum((tfmagNBand - target).^2, 'all')  / numel(tfmagNBand);
end