%% IIR filter loss function

function [loss, tfmagNBand, tfmag] = IIRFilterLoss(output, target, controlparameters)
    
    [z, p, k] = CreateIIRFromNNOutput(output, controlparameters.filterOrder);

    [tfmag, ~] = CreateIIRFilter(z, p, k, controlparameters.nfft, controlparameters.fs);

    tfmagNBand = CreateNBandMagnitude(tfmag, controlparameters.fidx);

    tfmagNBand = max(-128, min(128, tfmagNBand));
    
    %DC = max(-128, min(128, tfmag(1,:)));
    %target = target(2:end,:);
    %loss = 0.5 * mean((tfmagNBand - target).^2, 'all') + 0.5 * mean((tfmag(1,:) - DC).^2, 'all');
    %loss = mean([(tfmagNBand - target).^2; DC .^ 2], 'all');
    %loss = mean((tfmagNBand - target).^2, 'all') + mean(DC .^ 2, 'all');
    loss = mean((tfmagNBand - target).^2, 'all');
end