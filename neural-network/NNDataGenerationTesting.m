%close all
%clear all

z = [-0.9, -0.8];
p = [0.9, -0.9];
k = [0.1, 5];

fs = 48e3;
nfft = 8192;
[tfmag, fvec] = CreateIIRFilter(z, p, k, nfft, fs);

figure
semilogx(fvec, tfmag)
grid on
xlim([20 20e3])

%% Distributions
numObservations = 20e3;

weight = 25;
geometry = RandomGeometryWedge_Run2(numObservations);

wedgeIndex = geometry.wedgeIndex;
bendingAngle = geometry.bendingAngle;
minAngle = geometry.minAngle;
wedgeLength = geometry.wedgeLength;
radiusS = geometry.rS;
radiusR = geometry.rR;
zS = geometry.zS;
zR = geometry.zR;
zA = geometry.zA;

[~, phii] = CalculateApex(radiusS', radiusR', zS', zR', wedgeLength', true);
% zA = zA(:,3);

%% Plots

close all

zAProp = zA ./ wedgeLength;
figure
histogram(zAProp)
title('zA')

figure
histogram(zA(wedgeLength > 10), 'BinWidth', 1)
title('zA')

figure
histogram(zA, 'BinWidth', 1)
title('zA')

figure
histogram(phii, 'BinWidth', 1)
title('phii')

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

figure
histogram(wedgeLength, 'BinWidth', 1)
title('Wedge Length')

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

%%
figure
histogram(zS)
title('zS')

figure
histogram(zR)
title('zR')

test = [zS; zS];
figure
histogram([zS; zR])
title('zS and zR')

dZ = abs(zS - zR);
figure
histogram(dZ)
title('dZ')

%%
figure
histogram(zAProp(dZ < 0.5 & wedgeLength > 3))

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