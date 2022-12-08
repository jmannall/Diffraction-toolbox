%% Create output of IIR filters from ZPK parameters. Traceable by dlgradient

function [tfmag, fvec, tfcomplex] = CreateIIRFilter(z, p, k, nfft, fs)
    
    %% Calculate tfmag
    numObservations = size(z,2);

    numIIRFilters = size(z,1);
    
    [b, a] = IIRFilterCoefficients(z, p, k, numIIRFilters, numObservations);
    
    [tfmag, fvec, tfcomplex] = CalculateFilterResponse(b, a, nfft, fs);
end