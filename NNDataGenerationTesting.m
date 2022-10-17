clear all
close all

fs = 96e3;
nfft = 8192;

wedgeLength = 10;
wedgeIndex = 270;
thetaS = 20;
thetaR = 201;
radiusS = 1e6;
radiusR = 1;
[zS, zR] = deal(wedgeLength / 2);

controlparameters = struct('fs', fs, 'nfft', nfft, 'Rstart', radiusS);

epsilon = 1e-10;
createPlot = false;
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot);
[irDiffRef, tfmagDiffRef, tvecDiffRef, fvecDiffRef, tfcomplexDiffRef] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 + epsilon, radiusS, radiusR, zS, zR, controlparameters, createPlot);
[irDirRef, tfmagDirRef, tvecDirRef, fvecDirRef, tfcomplexDirRef] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 - epsilon, radiusS, radiusR, zS, zR, controlparameters, createPlot);

DiffRef = tfmagDiffRef.diff1;
DirRef = tfmagDirRef.direct;
shift = DirRef - DiffRef;

scaledResponse = [tfmag.diff1] + shift;
truth = 10;
i = min(1, (thetaR - thetaS - 180) / (wedgeIndex - thetaS - 180 - truth));

finalResponse = (1 - i) * scaledResponse + i * [tfmag.diff1];

figure
semilogx(fvecDiffRef, tfmagDiffRef.diff1, '--')
hold on
semilogx(fvecDirRef, tfmagDirRef.direct, '--')
semilogx(fvec, tfmag.diff1)
semilogx(fvec, scaledResponse)
semilogx(fvec, finalResponse)
legend('DiffRef', 'DirRef', 'Original', 'Scaled', 'Final', 'Location','southwest')
xlim([20 20e3])
ylim([-60 5])
