close all

wedge = 300;
bendingAngle = 5.5:359.5;
minAngle = 45.01;

geometry = Geometry(wedge, bendingAngle, minAngle);

disp('Run')
[result, geometry] = SingleWedgeArray(geometry, 20,1,1,10,10,48000);

tfmag = [result.tfmag];
fvec = result(1).fvec;
bendingAngle = geometry.bendingAngle';

%% Plot

numAngles = length(bendingAngle);
skip = 10;
extra = rem(numAngles, skip);

idxMin = floor(extra / 2) + 1;
idxMax = numAngles - ceil(extra / 2);
minbA = bendingAngle(idxMin);
maxbA = bendingAngle(idxMax);

logSkip = fvec(2) - fvec(1);

logIdx = round(logspace(log10(1), log10(4096), 1e5));
logFvec = fvec(logIdx);
logTfmag = tfmag(logIdx,:);

limits = [-70 0];

close all
figure
A = axes;
image(logTfmag,'CDataMapping','scaled')
A.CLim = limits;
colorbar
set(A,'YScale','log')
title('BTM target')
xticks(idxMin:skip:idxMax)
xticklabels(minbA:skip:maxbA)
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvec(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

%% NN
disp('Neural net')
load("NeuralNetwork_014.mat", 'net')

numBiquads = 2;
numFreq = 4096;

const = ones(1, numAngles);
%output = forward(net, dlarray([const * wedge; bendingAngle; const * minAngle; abs(bendingAngle - 180)], "CB"));
output = predict(net, dlarray([const * wedge; bendingAngle; const * minAngle], "CB"));

[zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
k = k / 100;
[tfmagNN, fvecNN] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

tfmagNN = extractdata(tfmagNN)';

logFvecNN = fvecNN(logIdx);
logTfmagNN = tfmagNN(logIdx,:);

figure
A = axes;
image(logTfmagNN,'CDataMapping','scaled')
A.CLim = limits;
colorbar
set(A,'YScale','log')
title('Neural network prediction')
xticks(idxMin:skip:idxMax)
xticklabels(minbA:skip:maxbA)
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvecNN(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])

error = abs(logTfmag - logTfmagNN);

figure
A = axes;
image(error,'CDataMapping','scaled')
A.CLim = [0 40];
colorbar
set(A,'YScale','log')
title('Neural network prediction error')
xticks(idxMin:skip:idxMax)
xticklabels(minbA:skip:maxbA)
%idxfvec = logspace(log10(1),log10(512),20);
%yticks(idxfvec)
%freq = round(fvec(round(idxfvec)));
yticklabels(round(logFvecNN(round(linspace(1, 1e5, 6)))))
ylim([0 1e5])