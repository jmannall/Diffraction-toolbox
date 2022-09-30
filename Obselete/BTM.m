%% Generate the BTM frequency response of an infinite wedge.

function [tfmag, fvec] = BTMInf(w, bA, mA, fs)
    
    % Inputs
    wedgeLength = 20;
    wedgeIndex = w;
    [radiusS, radiusR] = deal(1);
    [zS, zR] = deal(wedgeLength / 2);
    thetaS = mA;
    thetaR = mA + bA;

    % Create magnitude response and frequency vector
    [~, tfmag, ~, fvec, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, fs);
end