function H = SingleH(z, theta, f, gParameters, controlparameters)

    const = (2 * controlparameters.c) / (pi ^ 2 * gParameters.d * sin(gParameters.phii) ^ 2);

    %fc = const * CalculateNv(gParameters.v, theta) ^ 2;
    scale = CalculateInterpolationFactor(gParameters);
    fc = (controlparameters.c * (gParameters.v * sin(gParameters.v * pi)) ^ 2) / (2 * pi ^ 2 * gParameters.d * sin(gParameters.phii) ^ 2) * scale;
    dCorner = sqrt(gParameters.rS ^ 2 + (z - gParameters.zS) ^ 2) + sqrt(gParameters.rR ^ 2 + (z - gParameters.zR) ^ 2);
    t1 = dCorner / controlparameters.c;
    t = t1 - gParameters.t0;
    g = (2 / pi) * atan(pi * sqrt(2 * fc * t));

    fc = (1 / g ^ 2) * fc;

    H = g * SingleUDFA(f, fc, g);
end