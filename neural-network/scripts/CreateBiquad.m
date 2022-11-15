%% Create output of biquad filters from ZPK parameters. Traceable by dlgradient

function [tfmag, fvec, tfcomplex] = CreateBiquad(zR, zI, pR, pI, k, nfft, fs)
    
    %% Calculate tfmag
    numObservations = size(zR,2);

    numBiquads = size(zR,1);
    
    [b, a] = BiquadCoefficients(zR, zI, pR, pI, k, numBiquads, numObservations);
    
    [tfmag, fvec, tfcomplex] = CalculateFilterResponse(b, a, nfft, fs);
end

