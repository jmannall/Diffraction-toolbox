%% Biquad loss function

function loss = BiquadLoss(output, target, numBiquads, numFreq, fs, fidx)
    
    [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);

    [tfmag, ~] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs);
    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));
    loss = sum((tfmagNBand - target).^2, 'all')  / numel(tfmagNBand);
end