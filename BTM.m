function [tfmag, fvec] = BTM(w, bA, mA, fs)
    wedgeLength = 20;
    [radiusS, radiusR] = deal(1);
    [zS, zR] = deal(wedgeLength / 2);
    thetaS = mA;
    thetaR = mA + bA;
    wedgeIndex = w;
    [~, tfmag, ~, fvec, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, fs);
end