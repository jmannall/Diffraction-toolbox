%% IIR filter loss function

function [loss, tfmagNBand, tfmag] = IIRFilterLoss(output, target, numIIRFilters, nfft, fs, fidx)
    
    [z, p, k] = CreateIIRFromNNOutput(output, numIIRFilters);

    [tfmag, ~] = CreateIIRFilter(z, p, k, nfft, fs);
    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));

    loss = sum((tfmagNBand - target).^2, 'all')  / numel(tfmagNBand);
end