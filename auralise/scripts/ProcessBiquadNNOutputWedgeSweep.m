function [tfmag, tfcomplex, b, a] = ProcessBiquadNNOutputWedgeSweep(NN, wedgeIndex, bendingAngle, minAngle, numBiquads, pathLength, c, numFreq, fs)
    
    load(['data\NeuralNetwork_', NN, '.mat'], 'net')
    
    numPositions = length(bendingAngle);
    
    const = ones(1, numPositions);
    output = extractdata(predict(net, dlarray([const * wedgeIndex; bendingAngle; const * minAngle], "CB")));

    [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
    [tfmag, fvec, tfcomplex] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs);
    
    idx = bendingAngle < 180 - 2 * minAngle | (180 - minAngle < bendingAngle & bendingAngle < 180);
    tfcomplex = PhaseShift(tfcomplex, pathLength, fvec', c, idx);
    k(idx) = -k(idx);
    
    [b, a] = BiquadCoefficients(zR, zI, pR, pI, k, numBiquads, numPositions);
end