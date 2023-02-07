close all
%clear all

% dz = abs(trainingData(7,:) - trainingData(8,:));
% 
% figure
% histogram(dz)
% title('true dz')
% 
mA = trainingData(3,:);
% 
% figure
% histogram(mA)
% title('true mA')
% 
bA = trainingData(2,:);
% 
% figure
% histogram(bA)
% title('true bA')

s = [sin(mA); cos(mA)];
r = [sin(mA + bA); cos(mA + bA)];

figure
plot(s(1,:), s(2,:), 'o')
hold on
plot(r(1,:), r(2,:), 'x')

%% Distributions
numObservations = 20e3;

weight = 25;
geometry = RandomGeometryWedge(numObservations);

wedgeIndex = geometry.wedgeIndex;
bendingAngle = geometry.bendingAngle;
minAngle = geometry.minAngle;
wedgeLength = geometry.wedgeLength;
radiusS = geometry.rS;
radiusR = geometry.rR;
zS = geometry.zS;
zR = geometry.zR;
zA = geometry.zA;

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

[~,edges] = histcounts(radiusS);
figure
histogram(radiusS,edges)
title('radiusS')

[~,edges] = histcounts(radiusR);
figure
histogram(radiusR,edges)
title('radiusR')

figure
histogram(zS)
title('zS')

figure
histogram(zS)
title('zR')

figure
histogram([zS; zS])
title('zS and zR')

dZ = abs(zS - zR);
figure
histogram(dZ)
title('dZ')

zAProp = zA ./ wedgeLength;
figure
histogram(zAProp)
title('zA')

controlparameters = struct('fs', 96e3, 'nfft', 16384, 'difforder', 1, 'c', 344);

%% Figure
for i = 1:20
    L(i,1) = sqrt((radiusS(i) + radiusR(i)) ^ 2 + (zR(i) - zS(i)) ^ 2);
    [~, fvec, tfcomplex(i,:)] = SingleWedgeInterpolated(wedgeLength(i), wedgeIndex(i), minAngle(i), minAngle(i) + bendingAngle(i), radiusS(i), radiusR(i), zS(i), zR(i), controlparameters, false);
    tfcomplexDiffRef(i,:) = L(i) * tfcomplex(i,:);
    tfmagRef(i,:) = mag2db(abs(tfcomplexDiffRef(i,:)));
    controlparameters.Rstart = L(i);
    [~, fvec, tfcomplex(i,:)] = SingleWedgeInterpolated(wedgeLength(i), wedgeIndex(i), minAngle(i), minAngle(i) + bendingAngle(i), radiusS(i), radiusR(i), zS(i), zR(i), controlparameters, false);
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