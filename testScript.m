clear all
close all

fs = 96e3;
nfft = 8192;

wedgeLength = 10;
wedgeIndex = 270;
thetaS = 20;
thetaR = 260;
radiusS = 1;
radiusR = 1;
[zS, zR] = deal(wedgeLength / 2);

controlparameters = struct('fs', fs, 'nfft', nfft);
createPlot = false;

[~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot);


radiusS = 1e6;
radiusR = 2;
[tfmagInt, fvecInt, tfcomplexInt] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot);

% Add A and L to adjust the amplitudes. Check what need to train NN on and
% what need to do after to adjust the result (amplitude)

figure
semilogx(fvec, [tfmag.diff1])
hold on
%semilogx(fvec, mag2db(abs([tfcomplex.diff1])))
semilogx(fvecInt, tfmagInt)

PlotSpectogram([tfcomplex, tfcomplex], fvec, [0, 1], [-70 0], 'Test', false, false)
