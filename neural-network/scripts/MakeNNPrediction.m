function [tfmag, tfmagN] = MakeNNPrediction(net, input, pathLength, numFilters, fidx, controlparameters)

    output = extractdata(predict(net, input));
    [z, p, k] = CreateIIRFromNNOutput(output, numFilters);
    [tfmag, ~] = CreateIIRFilter(z, p, k, controlparameters.nfft, controlparameters.fs);
    tfmag = mag2db((1 ./ pathLength) .* db2mag(tfmag));    
    tfmagN = CreateNBandMagnitude(tfmag, fidx);
end