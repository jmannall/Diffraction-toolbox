function [tfmag, fvec, tfcomplex] = UTDSingleWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters)

    % Convert to radians
    thetaS = deg2rad(thetaS);
    thetaR = deg2rad(thetaR);
    wedgeIndex = deg2rad(wedgeIndex);
    phii = deg2rad(phii);

    c = controlparameters.c;
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    fvec = fs / nfft * (0:nfft / 2 - 1);
    k = 2 * pi * fvec/ c;
    n = wedgeIndex / pi;
    L = radiusR * radiusS / (radiusS + radiusR) * sin(phii) ^ 2;

    A = exp(-1i  * k * (radiusS + radiusR)) / sqrt(radiusS * radiusR * (radiusS + radiusR)); % From Pisha code. Not certain why this
    frontFactor = -A .* exp(-1i * (pi / 4)) ./ (2 * n * sqrt(2 * pi * k) * sin(phii));

    tfcomplex = frontFactor .* (EquationHalf(thetaR - thetaS, n, k, L) + EquationHalf(thetaR + thetaS, n, k, L));
    tfmag = mag2db(abs(tfcomplex));
end