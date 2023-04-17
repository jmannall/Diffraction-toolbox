function H = SingleH(z, theta, f, gParameters, controlparameters) % 79

    const = (2 * controlparameters.c) / (pi ^ 2 * gParameters.d * sin(gParameters.phii) ^ 2);
    fcOld = const * CalculateNv(gParameters.v, theta) ^ 2;
    scale = CalculateInterpolationFactor(gParameters); % 21
    fc = (controlparameters.c * (gParameters.v * sin(gParameters.v * pi)) ^ 2) / (2 * pi ^ 2 * gParameters.d * sin(gParameters.phii) ^ 2) * scale; % 11
    dCorner = sqrt(gParameters.rS ^ 2 + (z - gParameters.zS) ^ 2) + sqrt(gParameters.rR ^ 2 + (z - gParameters.zR) ^ 2); % 9 everytime + 2 once
    t1 = dCorner / controlparameters.c; % 1
    t = t1 - gParameters.t0; % 1
    g = (2 / pi) * atan(pi * sqrt(2 * fc * t)); % 7

    fc = (1 / g ^ 2) * fc; % 3

    H = g * SingleUDFA(f, fc, g); % 1 + n -> 24
end