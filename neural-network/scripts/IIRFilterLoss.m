%% IIR filter loss function

function [loss, tfmagNBand, tfmag] = IIRFilterLoss(output, target, numIIRFilters, nfft, fs, fidx)
    
    [z, p, k] = CreateIIRFromNNOutput(output, numIIRFilters);

    [tfmag, ~] = CreateIIRFilter(z, p, k, nfft, fs);

    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));
    
    %DC = target(1,:);
    %target = target(2:end,:);
    %loss = 0.5 * mean((tfmagNBand - target).^2, 'all') + 0.5 * mean((tfmag(1,:) - DC).^2, 'all');
    loss = mean((tfmagNBand - target).^2, 'all');
end