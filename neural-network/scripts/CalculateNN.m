function [tfcomplex, b, a] = CalculateNN(output, tfcomplex, validPath, pathLength, nfft, fs, biquad)
    
    numFilters = 2;
    numReceivers = length(pathLength);
    if biquad
        [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numFilters);
        [b, a] = BiquadCoefficients(zR, zI, pR, pI, k, numFilters, numReceivers);
        [~, ~, tfcomplexNN] = CreateBiquad(zR, zI, pR, pI, k, nfft, fs);
    else
        [z, p, k] = CreateIIRFromNNOutput(output, numFilters);
        [b, a] = IIRFilterCoefficients(z, p, k, numFilters, numReceivers);
        [~, ~, tfcomplexNN] = CreateIIRFilter(z, p, k, nfft, fs);
    end

    tfcomplex = tfcomplex(1:nfft / 2,:);
    tfcomplex(:,validPath) = extractdata(tfcomplexNN(:,validPath) ./ pathLength(validPath)');
end
