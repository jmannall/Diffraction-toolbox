close all
clear all

rng(1)
fs = 96e3;
nfft = 8192;
batchSize = 200;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'saveFiles', true);

wedgeLength = 10;
radiusS = 0.02;
radiusR = 0.03;
thetaS = 10;
thetaR = 179.9999 + thetaS;
wedgeIndex = 320;
zS = 5;
zR = 5;

[ir, tfmag, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
dirRef = tfmag.direct;
dirIr = ir.direct;
input = [zeros(10, 1); 1; zeros(44, 1)];
windowLength = 11;
pathLength = (radiusS + radiusR) * ones(1, length(input) / windowLength);
validPath = ones(size(pathLength));
c = 344;
[output, ir] = DelayLine(input, pathLength, windowLength, validPath, c, fs);
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
numObservations = 20e3;

geometry = RandomGeometryWedge(numObservations);

wedgeIndex = geometry.wedgeIndex;
bendingAngle = geometry.bendingAngle;
minAngle = geometry.minAngle;
wedgeLength = geometry.wedgeLength;
radiusS = geometry.radiusS;
radiusR = geometry.radiusR;
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
histogram(zS)
title('zS')

figure
histogram(zR)
title('zR')

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