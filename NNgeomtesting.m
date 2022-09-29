% NNSingleWedgeArray
clear all
close all
% NNGeometry --  check the single wedge convex case - to include w < 180

wedgeLength = 20;
wedgeIndex = 270;
thetaS = 30;
thetaR = 260;
radiusS = 1;
radiusR = 1;
zS = wedgeLength / 2;
zR = wedgeLength / 2;
fs = 48000;

[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);
wedgeLength = 50;
zS = wedgeLength / 2;
zR = wedgeLength / 2;
[ir, tfmag2, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);

f2 = figure(2);
movegui(f2,'north');
semilogx(fvec, tfmag)
hold on
semilogx(fvec, tfmag2)
xlim([20, 20000])
ylim([-40 0])

%% Create wedge array
close all
epsilon = 1e-3;
gParameters = struct('wedgeLength', [0.01, 50], 'wedgeIndex', [180 + epsilon, 360 - epsilon], ...
    'thetaS', epsilon, 'thetaR', epsilon, 'radiusS', [1, 1], 'radiusR', [1, 1], 'zS', 10, 'zR', 10, 'epsilon', epsilon);
numObservations = 10000;
fs = 48000;

[result, NNinput] = NNWedgeArray(gParameters, numObservations, fs);

findex = find([result(1).fvec] > 20000);