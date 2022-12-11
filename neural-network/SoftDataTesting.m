close all
clear all

fs = 96e3;
nfft = 16384;
c = 344;
testSize = 20e3;

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
[inputData, targetData, ~, fvec, fidx] = CreateBtmTrainingData(testSize, controlparameters, 1);
targetData = extractdata(targetData);

meanTargetData = mean(targetData, 2);

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
fv = 0.9995;
thresholdVar = fv * totalVar;
idx = cumVar < thresholdVar; % m = sum(idx)
numPos = sum(idx);
P = p(:,idx);

b = P' * (targetData - meanTargetData); % m first dimension. k second dimension
b = max(-3 * sqrt(llambda(idx)), min(3 * sqrt(llambda(idx)), b));

numFreq = length(fvec);
components = zeros(numFreq, testSize, numPos);
for i = 1:numPos
    components(:,:,i) = P(:,i) * b(i,:);
end

figure
semilogx(fvec, P(:,1))
hold on
semilogx(fvec, P(:,2))
semilogx(fvec, P(:,3))
legend('1', '2', '3')

softTargetData = meanTargetData + P * b;

n = 3;
m = 4;
for j = 1:5
    figure
    tiledlayout(n, m)
    for i = (j - 1) * n * m + 1:j * n * m
        nexttile
        semilogx(fvec, meanTargetData)
        hold on
        semilogx(fvec, targetData(:,i))
        semilogx(fvec, softTargetData(:,i))
        semilogx(fvec, squeeze(components(:,i,:)), '--')
        legend('mean', 'hard target', 'soft target', 'components', 'Location','best')
        ylim([-50 20])
        xlim([20 20e3])
    end
end

% figure
% semilogx(fvec, targetData)
% ylim([-70 0])
% xlim([20 20e3])
% title('Target data')
% 
% figure
% semilogx(fvec, softTargetData)
% ylim([-70 0])
% xlim([20 20e3])
% title('Soft target data')

%% IIR filter

lpFc = 1000;
hsFc = 20e3;
G = 0;
k = 0.5;
endTarget = 20e3;
fendidx = find(fvec > endTarget, 1);

[b, a] = IIRFilterParameterCoefficients(lpFc, hsFc, G, k, fs);

[tfmag, ~, ~] = CalculateFilterResponse(b, a, nfft, fs);
tfmagRef = CreateNBandMagnitude(tfmag, fidx);
tfmagRef = extractdata(tfmagRef);
tfmagRef = P(:,2);

scale = 46;
tfmagTarget = scale * tfmagRef;
lpFc = fvec(find(tfmagTarget - tfmagTarget(1) < -3, 1));
k = 10 ^ (P(1,2) * scale /  20);

[b, a] = LowPassCoefficients(lpFc, fs, k);
[tfmag, ~, ~] = CalculateFilterResponse(b, a, nfft, fs);
tfmag = CreateNBandMagnitude(tfmag, fidx);
tfmag = extractdata(tfmag);

G = tfmagTarget(fendidx) - tfmag(fendidx);
hsFc = (lpFc + endTarget) / 1.4;
[b, a] = IIRFilterParameterCoefficients(lpFc, hsFc, G, k, fs);

[tfmag, ~, ~] = CalculateFilterResponse(b, a, nfft, fs);
tfmagPredict = CreateNBandMagnitude(tfmag, fidx);
tfmagPredict = extractdata(tfmagPredict);

figure
semilogx(fvec, tfmagRef)
hold on
semilogx(fvec, tfmagTarget)
semilogx(fvec, tfmagPredict)
legend('Original', 'Target', 'Prediction')

%% input data

const = ones(1, testSize);
angleAxis = [0 * const; rad2deg(inputData(3,:)); rad2deg(sum(inputData(2:3,:))); rad2deg(inputData(1,:)); 360 * const];
distanceAxis = [-50 * const; -inputData(5,:); 0 * const; inputData(6,:); 50 * const];
heightAxis = [-50 * const; inputData(7,:); inputData(4,:); inputData(8,:); 50 * const];

figure
plot3(angleAxis(2:4,1:10), distanceAxis(2:4,1:10), heightAxis(2:4,1:10))
xlabel('angle')
ylabel('distance')
zlabel('height')
ax = gca;
