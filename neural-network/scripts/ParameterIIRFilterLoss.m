%% IIR filter loss function

function [loss, tfmagNBand] = ParameterIIRFilterLoss(output, target, nfft, fs, fidx)
    
    [lpFc, hsFc, G, k] = CreateFilterParametersFromNNOutput(output);

    [b, a] = deal(zeros(2, 2, length(lpFc)));
    [b(:,1,:), a(:,1,:)] = LowPassCoefficients(lpFc, fs, k);
    [b(:,2,:), a(:,2,:)] = HighShelfCoefficients(hsFc, G, fs);

    [tfmag, ~, ~] = CalculateFilterResponse(b, a, nfft, fs);
    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));

    loss = sum((tfmagNBand - target).^2, 'all')  / numel(tfmagNBand);
end