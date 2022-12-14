%% IIR filter loss function

function [loss, tfmagNBand] = ParameterIIRFilterLoss(output, target, nfft, fs, fidx)
    
    [lpFc, hsFc, G, k] = CreateFilterParametersFromNNOutput(output);

    [b, a] = IIRFilterParameterCoefficients(lpFc, hsFc, G, k, fs);

    [tfmag, ~, ~] = CalculateFilterResponse(b, a, nfft, fs);
    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));

    loss = sum((tfmagNBand - target).^2, 'all')  / numel(tfmagNBand);
end