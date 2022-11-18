close all
clear all


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