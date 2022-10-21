close all
clear all

numInputs = 1;
numOutputs = 9;
hiddenLayerSize = 2;
numLayers = 10;
alpha = 0.2;

Nn = [8, hiddenLayerSize * ones(1, numLayers)];
Cgx = 1;
Ni = numInputs;
No = numOutputs;
Nl = numLayers;

C = sum((Nn(1:Nl) + Cgx + 1) .* Ni .* Nn(2:Nl + 1)) + (Nn(Nl + 1) + 1) .* Ni .* No;

%% UTD

fs = 96e3;
nfft = 4096;
c = 344;

f = fs / nfft * (0:nfft / 2 - 1);
k = 2 * pi * f/ c;

wedgeIndex = 320;
thetaS = 30;
thetaR =319;
radiusS = 1;
radiusR = 1;
phii = 90; % B0

controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c);

[tfmagUTD, fvecUTD, tfcomplexUTD] = UTDSingleWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters);

[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(20, wedgeIndex, thetaS, thetaR, radiusS, radiusR, 10, 10, controlparameters, false);

figure
semilogx(fvecUTD, tfmagUTD)
hold on
semilogx(fvec, tfmag.diff1)
legend('UTD', 'BTM')
xlim([20 20e3])

