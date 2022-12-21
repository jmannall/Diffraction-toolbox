function [tfmag, fvec, tfcomplex] = SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters, A, B) % 174

    % Convert to radians
    thetaS = deg2rad(thetaS);
    thetaR = deg2rad(thetaR);
    wedgeIndex = deg2rad(wedgeIndex);
    phii = deg2rad(phii);

    B0 = sin(phii); % 1
    l = radiusS / B0;
    m = radiusR / B0;
    c = controlparameters.c;
    if isfield(controlparameters, 'fvec')
        fvec = controlparameters.fvec;
    else
        fvec = [125 500 2000 11050];
    end
    k = 2 * pi * fvec/ c; % Precompute
    n = wedgeIndex / pi; % 1

    if nargin > 7
        L = A * sin(phii)^2;
        frontFactor = -exp(-1i * (pi / 4)) ./ (2 * n * sqrt(2 * pi * k));
    else
        L = l * m  * B0^2 / (l + m); % 5
        A = exp(-1i  * k * (l + m)) / sqrt(l * m * (l + m)); % 9 - From Pisha code. Not certain why this
        frontFactor = -A .* exp(-1i * (pi / 4)) ./ (2 * n * sqrt(2 * pi * k) * B0); % 12
        B = 1;
    end

    tfcomplex = frontFactor .* (EquationHalf(thetaR - thetaS, n, k, L, B) + EquationHalf(thetaR + thetaS, n, k, L, B)); % 4 + 2 * eqHalf -> 146
    tfmag = mag2db(abs(tfcomplex));
    % tfmag = 20 * log10(sqrt(real(tfcomplex) ^ 2 + imag(tfcomplex) ^ 2));
    % 8
end