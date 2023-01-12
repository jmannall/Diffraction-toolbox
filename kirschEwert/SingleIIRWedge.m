%% Create IIR filter diffraction model Kirsch and Ewert

function [tfmag, b, a, fvec] = SingleIIRWedge(wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters)

    fs = controlparameters.fs;
    nfft = controlparameters.nfft;
    c = controlparameters.c;
    T = 1 / fs;
    
    [zA, phii] = CalculateApex(rS, rR, zS, zR, 100, true);
    zA = zA(:,3);
    
    bendingAngle = abs(thetaR - thetaS);
    minAngle = min(thetaS, wedgeIndex - thetaR);
    bendingAngle = deg2rad(bendingAngle) - pi;  % Adjust bending angle
    minAngle = deg2rad(minAngle);
    wedgeIndex = deg2rad(wedgeIndex);
    phii = deg2rad(phii);
    
    dS = sqrt(rS ^ 2 + (zA - zS) ^ 2);
    dR = sqrt(rR ^ 2 + (zR - zA) ^ 2);
    d = (2 * dS * dR) / (dS + dR);
    maxBendingAngle = wedgeIndex - (minAngle + pi);
    mp = 1 - 0.75 * tanh(1 / (2 * bendingAngle)) * sqrt(tanh(2 * minAngle));
    mw = (1 - 0.75 * (bendingAngle / maxBendingAngle) * sqrt(sin(wedgeIndex / 2)))^(-1); % sin(-wedgeIndex / 2) in the paper. But this leads to complex results. This gives a closer match to UTD too so probably a typo?.
    
    fc = (c * real(mw) * mp) / (3 * pi * d * (1 - cos(bendingAngle)) * sin(phii) ^ 2);
    %% LPF
    f0 = 1.11 * fc * (15.6 / fc) ^ 0.141;
    [bLpf, aLpf] = EDlpf(f0, T);
    
    %% Hsh
    actual = (log2(fs / 2) - log2(f0)) * 6;
    target = (log2(fs / 2) - log2(fc)) * 3;
    
    fsh =   209 * fc * (15.6 / fc) ^ 0.827;
    G = max(0, actual - target);
    [bHsh, aHsh] = EDhsh(fsh, G, T);
    
    pathLength = dS + dR;
    bLpf = 1 / pathLength * bLpf;
    b = [bLpf bHsh];
    a = [aLpf aHsh];
    [tfmag, fvec] = CalculateFilterResponse(b, a, nfft, fs);
end

