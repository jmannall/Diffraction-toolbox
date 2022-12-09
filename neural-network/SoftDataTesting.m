close all
clear all

fs = 96e3;
nfft = 16384;
c = 344;
testSize = 20e3;

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
[inputData, targetData, ~, fvec] = CreateBtmTrainingData(testSize, controlparameters, 1);
targetData = extractdata(targetData);

meanTargetData = mean(targetData, 2);
softTargetData = (0.9 * targetData + 0.1 * meanTargetData);

const = ones(size(targetData));

test = (targetData(:,1) - meanTargetData) * (targetData(:,1) - meanTargetData)';
test2 = (targetData(:,2) - meanTargetData) * (targetData(:,2) - meanTargetData)';
test3 = test + test2;
test4 = (targetData(:,1:2) - meanTargetData) * (targetData(:,1:2) - meanTargetData)';
S = 1 / (testSize - 1) * (targetData - meanTargetData) * (targetData - meanTargetData)';

C = cov(targetData');
[p, ~] = eig(C); % Column vectors
llambda = eig(C);

p = fliplr(p);
llambda = flipud(llambda);

totalVar = sum(llambda);

cumVar = cumsum(llambda);
fv = 0.99;
thresholdVar = fv * totalVar;
idx = cumVar < thresholdVar; % m = sum(idx)
P = p(:,idx);

b = P' * (targetData - meanTargetData); % m first dimension. k second dimension
b = max(-3 * sqrt(llambda(idx)), min(3 * sqrt(llambda(idx)), b));

softTargetData = meanTargetData + P * b;

n = 4;
m = 3;
tiledlayout(n, m)

for i = 1:n * m
    nexttile
    semilogx(fvec, meanTargetData)
    hold on
    semilogx(fvec, targetData(:,i))
    semilogx(fvec, softTargetData(:,i))
    legend('mean', 'hard target', 'soft target')
end

figure
semilogx(fvec, targetData)
ylim([-70 0])
xlim([20 20e3])
title('Target data')

figure
semilogx(fvec, softTargetData)
ylim([-70 0])
xlim([20 20e3])
title('Soft target data')