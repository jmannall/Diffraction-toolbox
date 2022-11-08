close all
clear all

rng(1)
fs = 96e3;
nfft = 8192;
c = 344;
batchSize = 200;
controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'difforder', 1, 'saveFiles', true);

wedgeLength = 10;
radiusS = 2;
radiusR = 3;
thetaS = 10;
thetaR = 200;
wedgeIndex = 320;
zS = 5;
zR = 5;

receiverPath = [-4 2 1.6
    1 2 1.6
    8 -4 1.6];
source = [1 -1 1.6];

updateRate = 20;

wedgeIndex = 270;
[receiverPath, receiverDirection, receiverData, speed] = CreatePath(receiverPath, updateRate, 'scene_1', wedgeIndex, source);

receivers = receiverPath;

radiusS = norm(source(1:2));
thetaS = wedgeIndex + sign(source(:,2)) .* acosd(source(:,1) ./ radiusS) - 180;
zS = source(3);

radiusR = vecnorm(receivers(:,1:2), 2, 2);
thetaR = wedgeIndex + sign(receivers(:,2)) .* acosd(receivers(:,1) ./ radiusR) - 180;
zR = receivers(:,3);

createPlot = false;
for i = 1:length(receivers)
    [tfmagI(:,i), fvecI, tfcomplexI(:,i)] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR(i), radiusS, radiusR(i), zS, zR(i), controlparameters, createPlot);
    [~, tfmag(:,i), ~, fvec, tfcomplex(:,i)] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR(i), radiusS, radiusR(i), zS, zR(i), controlparameters, createPlot);
end

tfdiffI = tfcomplexI;
tfdiff = [tfcomplex.complete];

PlotSpectrogram(tfdiff, fvec, 1:length(receivers), [-70 0], 'Btm scene_1', false, false, 'time')
PlotSpectrogram(tfdiffI, fvecI, 1:length(receivers), [-70 0], 'Interpolated scene_1', false, false, 'time')

figure
semilogx(fvec, [tfmag.diff1])
hold on
semilogx(fvecI, tfmagI)
xlim([20 20e3])
legend('BTM', 'BTMi')

[ir, tfmag, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
dirRef = tfmag.direct;
dirIr = ir.direct;
input = [1; zeros(5, 1)];
windowLength = 3;
pathLength = radiusS + radiusR;
validPath = 1;
c = 344;
[~, ir] = DelayLine(input, pathLength, windowLength, 1, c, fs);
tfmag = IrToTf(ir, nfft);
[trainingData, targetData, fvec, fNBand, fidx] = CreateBtmTrainingData(batchSize, controlparameters);

figure
semilogx(fc, extractdata(trainingData), '-')

x = 1;

% receiver = [-1 0 15];
% source = [0 0.5 -5];
% apex = [0 0 5];
% eVec = [0 0 1];
% 
% x = 0;
% y = 0;
% z = 1;
% 
% u = [x y z];
% u = u / norm(u);
% 
% bendingAngle = 270;
% theta = 180 - bendingAngle;
% c = cosd(theta);
% s = sind(theta);
% C = 1 - c;
% 
% d = [c -z * s y * s
%     z * s c -x * s
%     -y * s x * s c];
% Q = C * (u' * u) + d;
% 
% vReceiver = receiver * Q;
% 
% test1 = acosd(dot((vReceiver - source) / norm(vReceiver - source), eVec));
% test2 = acosd(dot((receiver - apex) / norm(apex - receiver), eVec));
%% Distributions
numObservations = 20e3;

geometry = RandomGeometryWedge(numObservations);

wedgeIndex = geometry.wedgeIndex;
bendingAngle = geometry.bendingAngle;
minAngle = geometry.minAngle;
wedgeLength = geometry.wedgeLength;
radiusS = geometry.rS;
radiusR = geometry.rR;
zS = geometry.zS;
zR = geometry.zR;
truezS = geometry.truezS;
truezR = geometry.truezR;
truezA = geometry.truezA;

%% Plots

close all

figure
histogram(wedgeIndex)
title('wedgeIndex')

figure
histogram(bendingAngle)
title('bendingAngle')

figure
histogram(minAngle)
title('minAngle')

[~,edges] = histcounts(log10(wedgeLength));
figure
histogram(wedgeLength,10.^edges)
set(gca, 'xscale','log')
title('wedgeLength')

[~,edges] = histcounts(log10(radiusS));
figure
histogram(radiusS,10.^edges)
set(gca, 'xscale','log')
title('radiusS')

[~,edges] = histcounts(log10(radiusR));
figure
histogram(radiusR,10.^edges)
set(gca, 'xscale','log')
title('radiusR')

figure
histogram(truezS)
title('zS')

figure
histogram(truezS)
title('zR')

figure
histogram([truezS; truezS])
title('zS and zR')

dZ = abs(truezS - truezR);
figure
histogram(dZ)
title('dZ')

if truezR < truezS
    sourcePart = radiusR ./ (radiusS + radiusR);
    zA = truezR + sourcePart .* dZ;
else
    sourcePart = radiusS ./ (radiusS + radiusR);
    zA = truezS + sourcePart .* dZ;
end

zAProp = truezA ./ wedgeLength;
figure
histogram(zAProp)
title('zA')

controlparameters = struct('fs', 96e3, 'nfft', 16384, 'difforder', 1);

%% Figure
for i = 1:20
    L(i,1) = sqrt((radiusS(i) + radiusR(i)) ^ 2 + (truezR(i) - truezS(i)) ^ 2);
    [~, fvec, tfcomplex(i,:)] = SingleWedgeInterpolated(wedgeLength(i), wedgeIndex(i), minAngle(i), minAngle(i) + bendingAngle(i), radiusS(i), radiusR(i), truezS(i), truezR(i), controlparameters, false);
    tfcomplexDiffRef(i,:) = L(i) * tfcomplex(i,:);
    tfmagRef(i,:) = mag2db(abs(tfcomplexDiffRef(i,:)));
    controlparameters.Rstart = L(i);
    [~, fvec, tfcomplex(i,:)] = SingleWedgeInterpolated(wedgeLength(i), wedgeIndex(i), minAngle(i), minAngle(i) + bendingAngle(i), radiusS(i), radiusR(i), truezS(i), truezR(i), controlparameters, false);
    tfcomplexDiff(i,:) = tfcomplex(i,:);
    tfmag(i,:) = mag2db(abs(tfcomplexDiff(i,:)));
    controlparameters = rmfield(controlparameters, 'Rstart');
end

figure
semilogx(fvec, tfmag)
grid on
hold on
semilogx(fvec, tfmagRef, '--')
xlim([20 20e3])
ylim([-100 10])