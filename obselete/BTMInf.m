%% Generate the BTM frequency response of an infinite wedge.

function [tfmag, fvec] = BTMInf(S)
    
    % Inputs
    wedgeLength = 20;
    wedgeIndex = S.w;
    [radiusS, radiusR] = deal(1);
    [zS, zR] = deal(wedgeLength / 2);
    thetaS = S.mA;
    thetaR = S.mA + S.bA;
    fs = S.fs;

    % Create magnitude response and frequency vector
    [~, tfmag, ~, fvec, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, fs);
    tfmag = [tfmag.diff];
end