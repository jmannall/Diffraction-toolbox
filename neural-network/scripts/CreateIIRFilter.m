%% Create output of IIR filters from ZPK parameters. Traceable by dlgradient

function [tfmag, fvec, tfcomplex] = CreateIIRFilter(z, p, k, nfft, fs)
    
    %% Calculate tfmag
    [numIIRFilters, numObservations] = size(z);
    
    [b, a] = IIRFilterCoefficients(z, p, k, numIIRFilters, numObservations);
    
    [tfmag, fvec, tfcomplex] = CalculateFilterResponse(b, a, nfft, fs);
end