function [tfmag, fvec, tfcomplex] = SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters, A, B) % 174

    % Convert to radians
    thetaS = deg2rad(thetaS);
    thetaR = deg2rad(thetaR);
    wedgeIndex = deg2rad(wedgeIndex);
    phii = deg2rad(phii);

    c = controlparameters.c;
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    fvec = fs / nfft * (0:nfft / 2 - 1);
    k = 2 * pi * fvec/ c; % Precompute
    n = wedgeIndex / pi; % 1

    if nargin > 7
        L = A * sin(phii)^2;
        frontFactor = -exp(-1i * (pi / 4)) ./ (2 * n * sqrt(2 * pi * k));
    else
        B0 = sin(phii); % 1
        L = radiusR * radiusS / (radiusS + radiusR) * B0^2; % 5
        A = exp(-1i  * k * (radiusS + radiusR)) / sqrt(radiusS * radiusR * (radiusS + radiusR)); % 9 - From Pisha code. Not certain why this
        frontFactor = -A .* exp(-1i * (pi / 4)) ./ (2 * n * sqrt(2 * pi * k) * B0); % 12
        B = 1;
    end

    tfcomplex = frontFactor .* (EquationHalf(thetaR - thetaS, n, k, L, B) + EquationHalf(thetaR + thetaS, n, k, L, B)); % 4 + 2 * eqHalf -> 146
    tfmag = mag2db(abs(tfcomplex));
end