close all

fs = 96e3;
nfft = 4096;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 3);

numX = 199;
numY = 199;
numReceivers = numX * numY;
x = linspace(-10, 10, numX)';
y = linspace(-10, 10, numY)';

count = 0;
receiver = zeros(numReceivers, 3);
for i = 1:numX
    for j = 1:numY
        count = count + 1;
        receiver(count,:) = [x(i) y(j) 5];
    end
end
numPillars = 13;
gratingWidth = 0.5;
pillarHeight = 10;
createPlot = true;

n = (numPillars - 1) / 2;
pillarCenters = 3 * gratingWidth * [-fliplr(1:n), 0, 1:n];
xDist = min(abs(receiver(:,1) - pillarCenters), [], 2) <= gratingWidth;
yDist = abs(receiver(:,2)) <= gratingWidth;
inPillar = xDist & yDist;
count = sum(inPillar);
toCalculate = numReceivers - count;
parfor i = 1:numReceivers
    if inPillar(i)
        test = receiver(i,:);
        [tfmag(:,i), fvec(:,i), tfcomplex(:,i)] = deal(zeros(3, 1));
    else
        [tfmag(:,i), fvec(:,i), tfcomplex(:,i)] = SingleDiffractionGrating(numPillars, gratingWidth, pillarHeight, receiver(i,:), controlparameters, createPlot);
    end
end
fvec = fvec(:,1);

tfcomplex = reshape(tfcomplex, 3, [], numX);
%% Figure
PlotSpectrogram(squeeze(tfcomplex(1,:,:))', x, y, [-70 0], [num2str(fvec(1)), ' Hz'], false, false, 'x')
PlotSpectrogram(squeeze(tfcomplex(2,:,:))', x, y, [-70 0], [num2str(fvec(2)), ' Hz'], false, false, 'x')
PlotSpectrogram(squeeze(tfcomplex(3,:,:))', x, y, [-70 0], [num2str(fvec(3)), ' Hz'], false, false, 'x')



