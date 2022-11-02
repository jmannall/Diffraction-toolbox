close all
clear all

fs = 96e3;
nfft = 4096;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 3);

numX = 50;
numY = 50;
numReceivers = numX * numY;
epsilon = 1e-4;
x = linspace(-1.5 + epsilon, 1.5 - epsilon, numX)';
y = linspace(-1.5 + epsilon, 1.5 - epsilon, numY)';

count = 0;
receiver = zeros(numReceivers, 3);
for i = 1:numX
    for j = 1:numY
        count = count + 1;
        receiver(count,:) = [x(i) y(j) 1.5];
    end
end
numPillars = 13;
gratingWidth = 0.3 / 3.5;
pillarWidth = 0.5 / 3.5;
pillarHeight = 3;
createPlot = false;

[result, geometry] = SingleDiffractionGratingArray(numPillars, gratingWidth, pillarWidth, pillarHeight, receiver, controlparameters, createPlot);

fvec = result(1).fvec;
tfcomplex = [result.tfcomplex];

tfcomplex = reshape([tfcomplex.complete], 3, [], numX);

PlotSpectrogram(squeeze(tfcomplex(1,:,:)), y, x, [-70 0], [num2str(fvec(1)), ' Hz'], false, false, 'x')
PlotSpectrogram(squeeze(tfcomplex(2,:,:)), y, x, [-70 0], [num2str(fvec(2)), ' Hz'], false, false, 'x')
PlotSpectrogram(squeeze(tfcomplex(3,:,:)), y, x, [-70 0], [num2str(fvec(3)), ' Hz'], false, false, 'x')



