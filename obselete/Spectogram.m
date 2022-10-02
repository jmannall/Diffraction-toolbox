clear all
close all

n = 798;
wedge = 350;
minAngle = 30.01;
bendingAngle = 4.99:359.99;
bendingAngle = linspace(4.99, wedge - minAngle - 0.01, n);

geometry = GeometryAll(wedge, bendingAngle, minAngle);

fs = 48000;
disp('Run')
[result, geometry] = SingleWedgeArray(geometry, 20,1,1,10,10,fs);

tfmag = [result.tfmag];
tfcomplex = [result.tfcomplex];
fvec = result(1).fvec;
bendingAngle = geometry.bendingAngle';

%% Plot

numAngles = length(bendingAngle);
skip = 20;
extra = rem(numAngles, skip);

idxMin = floor(extra / 2) + 1;
idxMax = numAngles - ceil(extra / 2);
minbA = bendingAngle(idxMin);
maxbA = bendingAngle(idxMax);

logSkip = fvec(2) - fvec(1);

logIdx = round(logspace(log10(1), log10(4096), 1e5));
logFvec = fvec(logIdx);
tfmagComplete = [tfmag.complete];
tfmagDirect = [tfmag.direct];
tfmagGeom = [tfmag.geom];
tfmagDiff = [tfmag.diff];

tfphaseDiff = angle([tfcomplex.diff]);

logTfmagComplete = tfmagComplete(logIdx,:);
logTfmagDirect = tfmagDirect(logIdx,:);
logTfmagGeom = tfmagGeom(logIdx,:);
logTfmagDiff = tfmagDiff(logIdx,:);
logTfphaseDiff = tfphaseDiff(logIdx,:);

magLimits = [-70 0];
phaseLimits = [-pi pi];

close all

figure
A = axes;
image(logTfmagComplete,'CDataMapping','scaled');
A.CLim = magLimits;
colorbar
set(A,'YScale','log')
title('BTM Complete')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvec(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

figure
A = axes;
image(logTfmagDirect,'CDataMapping','scaled');
A.CLim = magLimits;
colorbar
set(A,'YScale','log')
title('BTM Direct')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvec(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

figure
A = axes;
image(logTfmagGeom,'CDataMapping','scaled');
A.CLim = magLimits;
colorbar
set(A,'YScale','log')
title('BTM Specular')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvec(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

figure
A = axes;
image(logTfmagDiff,'CDataMapping','scaled');
A.CLim = magLimits;
colorbar
set(A,'YScale','log')
title('BTM Diffraction magnitude')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvec(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

figure
A = axes;
image(logTfphaseDiff,'CDataMapping','scaled');
%A.CLim = phaseLimits;
colorbar
set(A,'YScale','log')
title('BTM Diffraction phase')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvec(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

%% NN
disp('Neural net')
load("NeuralNetwork_015b.mat", 'net')

numBiquads = 2;
numFreq = 4096;

const = ones(1, numAngles);
output = predict(net, dlarray([const * wedge; bendingAngle; const * minAngle], "CB"));

[zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);

idx = bendingAngle < 180 - 2 * minAngle | (180 - minAngle < bendingAngle & bendingAngle < 180);
k(idx) = -k(idx);

[tfmagNN, fvecNN, tfcomplexNN] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs);

z = [extractdata(zR) + 1i*extractdata(zI)];
p = [extractdata(pR) + 1i*extractdata(pI)];
k = extractdata(k);

fs = 48000;
c = 343;
d = 2;
ts = 1/ fs;
delay = d / c;
sampleDelay = round(delay / ts);

mag = abs(tfcomplexNN);
phase = angle(tfcomplexNN);
phaseShift = -2 * pi * delay * fvec';
phase = phase + phaseShift;

tfcomplexNN = mag.*exp(phase*1i);
tfcomplexNN = extractdata(tfcomplexNN);

tfmag = 20*log10(abs(tfcomplexNN));
tfmagNN = extractdata(tfmagNN)';

logFvecNN = fvecNN(logIdx);
logTfmagNN = tfmagNN(logIdx,:);

tfphaseNN = angle(tfcomplexNN);
logTfphaseNN = tfphaseNN(logIdx,:);

tfcomplexNNAll = tfcomplexNN + [tfcomplex.direct] + [tfcomplex.geom];
tfmagNNAll = 20*log10(abs(tfcomplexNNAll));
logTfmagNNAll = tfmagNNAll(logIdx,:);

figure
A = axes;
image(logTfmagNN,'CDataMapping','scaled')
A.CLim = magLimits;
colorbar
set(A,'YScale','log')
title('Neural network prediction')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
yticklabels(round(logFvecNN(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

figure
A = axes;
image(logTfmagNNAll,'CDataMapping','scaled')
A.CLim = magLimits;
colorbar
set(A,'YScale','log')
title('Neural network prediction magnitude')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
yticklabels(round(logFvecNN(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

figure
A = axes;
image(logTfphaseNN,'CDataMapping','scaled')
%A.CLim = phaseLimits;
colorbar
set(A,'YScale','log')
title('Neural network prediction phase')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
yticklabels(round(logFvecNN(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

error = abs(logTfmagComplete - logTfmagNNAll);

figure
A = axes;
image(error,'CDataMapping','scaled')
A.CLim = [0 10];
colorbar
set(A,'YScale','log')
title('Neural network prediction error')
xticks(idxMin:skip:idxMax)
xticklabels(round(minbA:skip:maxbA))
yticklabels(round(logFvecNN(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])