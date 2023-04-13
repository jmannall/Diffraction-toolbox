function [tfmag, tfmagN, tfcomplex] = MakeNNPrediction(net, input, pathLength, numFilters, fidx, controlparameters)

    output = extractdata(predict(net, input));
    [z, p, k] = CreateIIRFromNNOutput(output, numFilters);
    [~, ~, tfcomplex] = CreateIIRFilter(z, p, k, controlparameters.nfft, controlparameters.fs);
    tfcomplex = (1 ./ pathLength) .* tfcomplex;
    tfmag = mag2db(abs(tfcomplex));    
    tfmagN = CreateNBandMagnitude(tfmag, fidx);
end