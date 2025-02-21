function [tfmag, tfmagN, tfcomplex] = MakeUDFA_NNPrediction(net, input, controlparameters)

    output = extractdata(predict(net, input));
    [z, p, k] = CreateIIRFromNNOutput(output, controlparameters.filterOrder);
    [~, ~, tfcomplex] = CreateIIRFilter(z, p, k, controlparameters.nfft, controlparameters.fs);
    tfmag = mag2db(abs(tfcomplex));    
    tfmagN = CreateNBandMagnitude(tfmag, controlparameters.fidx);
end