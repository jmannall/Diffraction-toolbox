function H = CalculateUDFATarget(f, gParameters, controlparameters) % 336

    z = gParameters.wL;
    [zA, phii] = CalculateApex(gParameters.rS, gParameters.rR, gParameters.zS, gParameters.zR, 0, true);
    gParameters.zA = zA(:,3);
    gParameters.phii = deg2rad(phii);

    gParameters.dS = sqrt(gParameters.rS ^ 2 + (gParameters.zA - gParameters.zS) ^ 2); % 5
    gParameters.dR = sqrt(gParameters.rR ^ 2 + (gParameters.zR - gParameters.zA) ^ 2); % 5
    gParameters.d = (2 * gParameters.dS * gParameters.dR) / (gParameters.dS + gParameters.dR); % 4
    gParameters.t0 = (gParameters.dS + gParameters.dR) / controlparameters.c; % 1

    gParameters.v = pi / gParameters.wI; % 1
    H = (SingleHpm(0, f, gParameters, controlparameters) + SingleHpm(z, f, gParameters, controlparameters)) / 4; % 2n + 2 -> 320
end