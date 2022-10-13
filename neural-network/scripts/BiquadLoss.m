%% Biquad loss function

function [loss, delloss] = BiquadLoss(output, target, numBiquads, numFreq, fidx)
    
    [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);

    [tfmag, ~] = CreateBiquad(zR, zI, pR, pI, k, numFreq);
    
    tfmagNBand = CreateNBandMagnitude(tfmag', fidx);

    loss = (tfmagNBand - 128 * target).^2 / (numel(tfmagNBand));
    delloss = dlgradient(sum(loss, "all"), output);
end