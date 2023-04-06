
close all
clear all
set(0, 'DefaultLineLineWidth', 1.5);

%% Constants

fs = 48e3;
nfft = 8192;
c = 344;
testSize = 20e3;
nBands = 8;

colorStore = colororder;

%% NN Parameters

%netName = 'iir-2057_0001-1-09-099-3-25.mat';
netName = ['Run2', filesep, 'IIR-6_24_0001.mat'];

%% Geometry input data

disp('Load validation data')

controlparameters = struct('fs', 2 * fs, 'nfft', 2 * nfft, 'difforder', 1, 'c', c, 'saveFiles', 3, 'noDirect', false);
[inputData, targetData, fvec] = CreateBtmTrainingData(testSize, controlparameters, 'ValidationData');
numFreq = length(fvec);

controlparameters.fs = fs;
controlparameters.nfft = nfft;
[~, tfmag.Default, ~, fvec] = DefaultBTM(controlparameters);

[~, fc, fidx, fLow, fHigh, BW] = CreateFrequencyNBands(tfmag.Default, fvec, nBands);

AWeighting = weightingFilter('A-weighting',fs);
aWeight = freqz(AWeighting, fc, fs);
aWeight = abs(aWeight)';

wI = inputData(1,:);
bA = inputData(2,:);
mA = inputData(3,:);
wL = inputData(4,:);
rS = inputData(5,:);
rR = inputData(6,:);
zS = inputData(7,:);
zR = inputData(8,:);

thetaS = rad2deg(mA);
thetaR = rad2deg(mA + bA);
wedgeIndex = rad2deg(wI);

[zA, phii] = CalculateApex(rS, rR, zS, zR, wL, false);
dZ = abs(zR - zS);
pathLength = sqrt((rS + rR) .^ 2 + dZ .^ 2);

%% BTM

disp('BTM')

% tfmag.BtmI = zeros(length(fvec), testSize);
% for i = 1:testSize
%     tfmag.BtmI(:,i) = SingleWedgeInterpolated(wL(i), wedgeIndex(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), controlparameters, false);
%     [~, tfmag.Btm(:,i)] = SingleWedge(wL(i), wedgeIndex(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), controlparameters, false);
% end
% tfmag.BtmI = mag2db((1 ./ pathLength) .* db2mag(tfmag.BtmI));
% tfmag.Btm = [tfmag.Btm.diff1];
% [tfmagN.BtmI, ~, ~] = CreateFrequencyNBands(tfmag.BtmI, fvec, nBands);
tfmagN.BtmI = mag2db((1 ./ pathLength) .* db2mag(targetData));

%% UTD LR

disp('UTD LR')


tf = zeros(4, testSize);
if isfield(controlparameters, 'fvec')
    controlparameters = rmfield(controlparameters, 'fvec');
end
for i = 1:testSize
    [tf(:,i), ~, ~] = SingleUTDWedgeInterpolated(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii(i), controlparameters);
end

validPath = true(size(pathLength));
input = [1; zeros(5, 1)];
windowLength = length(input);

fileName = 'UtdTestSetLR_irs.mat';
fileExists = exist(fileName, "file");
if fileExists == 2
    load("UtdTestSetLR_irs.mat", 'ir')
else
    [~, ir] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tf, c, fs, validPath);
    fileName = extractBefore(fileName, '.');
    save(fileName, 'ir')
end
tfmag.UtdILR = IrToTf(ir, nfft);
fvecUtd = fvec(2:end);
[tfmagN.UtdILR, ~, ~] = CreateFrequencyNBands(tfmag.UtdILR, fvecUtd, nBands);

%% UTD

disp('UTD')

controlparameters.fvec =  fvecUtd;
tfmag.UtdI = zeros(length(fvecUtd), testSize);
for i = 1:testSize
    tfmag.UtdI(:,i) = SingleUTDWedgeInterpolated(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii(i), controlparameters);
end
[tfmagN.UtdI, ~, ~] = CreateFrequencyNBands(tfmag.UtdI, fvecUtd, nBands);

%% NN

disp('NN')

loadDir = 'NNSaves';
load([loadDir, filesep, netName], 'net');

numFilters = 2;
X = dlarray(single(inputData), "CB");
[tfmag.NN, tfmagN.NN] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);


%% IIR filter Kirsch and Ewert

disp('IIR')

tfmag.IIR = zeros(length(fvec), testSize);
for i = 1:testSize
    tfmag.IIR(:,i) = SingleIIRWedge(wedgeIndex(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), controlparameters);
end
tfmagN.IIR = CreateFrequencyNBands(tfmag.IIR, fvec, nBands);

%% Plots constants

savePath = 'figures';
fontSize = 18;
x = 300;
gap = 10;
numSets = 5;
n = 3;
m = 4;

roundInput = round([rad2deg(inputData(1:3,:)); inputData(4:8,:)], 1);

%% Frequency Plots

disp('Frequency Plots')

close all

for j = 1:numSets
    figure('Position', [gap gap 6 * x 3 * x])
    t = tiledlayout(n ,m, 'TileSpacing', 'compact');
    title(t, 'Frequency responses for given input parameters')
    for i = 1:n * m
        idx = (j - 1) * n * m + i;
        nexttile
        semilogx(fc, tfmagN.BtmI(:,idx))
        hold on
        semilogx(fvec, tfmag.NN(:,idx))
        semilogx(fvecUtd, tfmag.UtdI(:,idx))
        semilogx(fvec, tfmag.UtdILR(:,idx))
        semilogx(fvec, tfmag.IIR(:,idx))
        grid on
        xlim([50 20e3])
        ylim([-70 0])
        l = legend('BTMi', 'NN-IIR', 'UTDi-LR', 'UTDi', 'IIR');
        title(l, 'Diffraction Model')
        ylabel('Magnitude (dB)')
        xlabel('Frequency (Hz)')
        title({['wL: ' num2str(roundInput(4, idx)), ' wI: ', num2str(roundInput(1, idx)), ' mA: ', num2str(roundInput(3, idx)), ' bA: ', num2str(roundInput(2, idx))] ...
            [' rS: ' num2str(roundInput(5, idx)), ' rR: ', num2str(roundInput(6, idx)), ' zS: ', num2str(roundInput(7, idx)), ' zR: ', num2str(roundInput(8, idx))]})
    end
    %saveas(gcf, [savePath, filesep, 'AccuracyExample_', num2str(j), '.svg'], 'svg')
end

%% Data

close all

lossInf = CalculateLoss(tfmagN, tfmagN.BtmI, wL > 5);

figure
semilogx(fc, lossInf.f.NN)
hold on
semilogx(fc, lossInf.f.UtdILR)
semilogx(fc, lossInf.f.UtdI)
semilogx(fc, lossInf.f.IIR)
grid on
legend('NN', 'UtdILR', 'UtdI', 'IIR')
xlim([20 20e3])
ylim([0 12])

lossN = CalculateLoss(tfmagN, tfmagN.BtmI);

figure
semilogx(fc, lossN.f.NN)
hold on
semilogx(fc, lossN.f.UtdILR)
semilogx(fc, lossN.f.UtdI)
semilogx(fc, lossN.f.IIR)
grid on
legend('NN', 'UtdILR', 'UtdI', 'IIR')
xlim([20 20e3])
ylim([0 12])

%% Max loss case

[maxLoss, idx] = max(lossN.i.NN);
figure('Position', [gap gap 2 * x x])
semilogx(fc, tfmagN.BtmI(:,idx))
hold on
semilogx(fvec, tfmag.NN(:,idx))
semilogx(fvec, tfmag.UtdILR(:,idx))
semilogx(fvecUtd, tfmag.UtdI(:,idx))
semilogx(fvec, tfmag.IIR(:,idx))
grid on
xlim([50 20e3])
ylim([-70 0])
l = legend('BTM', 'NN-IIR', 'UTD-LR', 'UTD', 'IIR');
title(l, 'Diffraction Model')
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
title({['MAX wL: ' num2str(roundInput(4, idx)), ' wI: ', num2str(roundInput(1, idx)), ' mA: ', num2str(roundInput(3, idx)), ' bA: ', num2str(roundInput(2, idx))] ...
    [' rS: ' num2str(roundInput(5, idx)), ' rR: ', num2str(roundInput(6, idx)), ' zS: ', num2str(roundInput(7, idx)), ' zR: ', num2str(roundInput(8, idx))]})

%% Histograms

limit = [0 12e3];
figure('Position', [gap gap 4 * x 3 * x])
t = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t, 'Losses compared to BTM')
nexttile
histogram(lossN.i.UtdI)
title('UTD losses')
xlabel('L1 norm loss')
ylabel('Number of cases')
ylim(limit)
grid on

[~,edges] = histcounts(lossN.i.UtdI);

nexttile
histogram(lossN.i.UtdILR, edges)
title('UTD-LR losses')
xlabel('L1 norm loss')
ylabel('Number of cases')
ylim(limit)
grid on

nexttile
histogram(lossN.i.NN, edges)
title('NN-IIR losses')
xlabel('L1 norm loss')
ylabel('Number of cases')
ylim(limit)
grid on

nexttile
histogram(lossN.i.IIR, edges)
title('IIR losses')
xlabel('L1 norm loss')
ylabel('Number of cases')
ylim(limit)
grid on

%saveas(gcf, [savePath, filesep, 'AccuracyHistograms.svg'], 'svg')


%% Shaded plots

close all

percentiles = 10:10:50;

figure
PlotVarDependentLoss(lossN.i.NN, wL, 1, percentiles, 'wL', 1)
PlotVarDependentLoss(lossN.i.UtdI, wL, 1, percentiles, 'wL', 2)


figure
PlotVarDependentLoss(lossInf.i.NN, phii(wL > 5), 1, percentiles, 'Phii', 1)
PlotVarDependentLoss(lossInf.i.UtdI, phii(wL > 5), 1, percentiles, 'Phii', 2)

figure
PlotVarDependentLoss(lossInf.i.NN, rad2deg(mA(wL > 5)), 1, percentiles, 'mA', 1)
PlotVarDependentLoss(lossInf.i.UtdI, rad2deg(mA(wL > 5)), 1, percentiles, 'mA', 2)

figure
PlotVarDependentLoss(lossInf.i.NN, rad2deg(bA(wL > 5)), 1, percentiles, 'bA', 1)
PlotVarDependentLoss(lossInf.i.UtdI, rad2deg(bA(wL > 5)), 1, percentiles, 'bA', 2)

figure
PlotVarDependentLoss(lossInf.i.NN, wedgeIndex(wL > 5), 1, percentiles, 'wI', 1)
PlotVarDependentLoss(lossInf.i.UtdI, wedgeIndex(wL > 5), 1, percentiles, 'wI', 2)

figure
PlotVarDependentLoss(lossInf.i.NN, rS(wL > 5) + rR(wL > 5), 1, percentiles, 'L', 1)
PlotVarDependentLoss(lossInf.i.UtdI, rS(wL > 5) + rR(wL > 5), 1, percentiles, 'L', 2)

figure('Color','w');
subplot(3,1,1);
plot(X,Y,'LineWidth',1.5);
title('plot (True value y=f(x))');
ylim([-1 5]);

subplot(3,1,2);
hold on
plot(X,Y,'LineWidth',1.5);
plot_distribution(X,Y_noisy);
hold off
title('plot\_distribution');
ylim([-1 5]);

subplot(3,1,3);
hold on
plot(X,Y,'LineWidth',1.5);
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
hold off
title('plot\_distribution\_prctile');
ylim([-1 5]);





idx = wL > 0;
test = sum(idx);

percentiles = 0:0.1:100;

percentilesNN1Sq = prctile(iSqLossesNN1(idx), percentiles);
percentilesNN2Sq = prctile(iSqLossesNN2(idx), percentiles);

percentilesUtd = prctile(iLossesUtd(idx), percentiles);
percentilesUtdLR = prctile(iLossesUtdLR(idx), percentiles);
percentilesNN1 = prctile(iLossesNN1(idx), percentiles);
percentilesNN2 = prctile(iLossesNN2(idx), percentiles);
percentilesIIR = prctile(iLossesIIR(idx), percentiles);

color = colorStore([1:3, 5],:);

figure
plot(percentiles, percentilesNN1)
hold on
plot(percentiles, percentilesUtdLR, '--')
plot(percentiles, percentilesUtd, '-.')
plot(percentiles, percentilesNN2)
%plot(percentiles, percentilesIIR, ':')
grid on
colororder(color)
l = legend('NN-IIR', 'UTD-LR', 'UTD', 'NN-IIR new', 'Location','northwest');
title(l, 'Diffraction Model')
xlabel('Percentile')
ylabel('Mean absolute error (dBA)')
ylim([0 30])
%title('Value at increasing percentiles of the L1 norm loss (A-weighted)')

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.2,0.2,0.35,0.4])
box(ax,'on')
plot(percentiles, percentilesNN1, 'parent',ax)
hold on
plot(percentiles, percentilesUtdLR,  '--', 'parent',ax)
plot(percentiles, percentilesUtd, ' -.', 'parent',ax)
plot(percentiles, percentilesNN2, 'parent',ax)
%plot(percentiles, percentilesIIR, ':', 'parent',ax)
grid on
set(ax,'xlim',[0,40],'ylim',[0,1])

%saveas(gcf, [savePath, filesep, 'PercentileLoss'], 'epsc')
%saveas(gcf, [savePath, filesep, 'PercentileLoss'], 'svg')

figure
plot(percentiles, percentilesNN1)
hold on
plot(percentiles, percentilesUtdLR)
plot(percentiles, percentilesUtd)
plot(percentiles, percentilesNN2)
%plot(percentiles, percentilesIIR)
grid on
legend('NN-IIR', 'UTD-LR', 'UTD', 'NN-IIR new' , 'Location','northwest')
xlabel('Percentile')
ylabel('Mean absolute error (dBA)')
title('Value at increasing percentiles of the L1 norm loss (A-weighted)')
xlim([0 40])

%saveas(gcf, [savePath, filesep, 'PercentileLoss_Cropped.svg'], 'svg')

%% Score

%close all

score = sum(iLossesNN2 < iLossesUtdLR) / testSize * 100;
comparison = iLossesUtdLR - iLossesNN2;

figure
bar(comparison)
grid on
title('Improvement in NN-IIR over UTD-LR')
ylabel('Mean absolute error improvement (dBA)')
xlabel('Number of cases')
ylim([-10 35])
xticks(0:2e3:2e4)
xticklabels(0:2e3:2e4)

%saveas(gcf, [savePath, filesep, 'Comparison.svg'], 'svg')

percentilesComparison = prctile(comparison, percentiles);

figure
plot(percentiles, percentilesComparison)
grid on
xlabel('Mean absolute error improvement (dBA)')
ylabel('Percentile')

figure
histogram(comparison)
grid on
%title('Improvement in NN-IIR over UTD-LR')
xlabel('Mean absolute error improvement (dBA)')
ylabel('Number of cases')
xlim([-10 35])

%saveas(gcf, [savePath, filesep, 'ComparisonHistogram'], 'epsc')
%saveas(gcf, [savePath, filesep, 'ComparisonHistogram'], 'svg')

%% Function of wedgeLength

%close all

[x, idx] = sort(inputData(4,:));
UtdLR = iLossesUtdLR(idx);
NN = iLossesNN2(idx);
y = [NN; UtdLR];

wL = inputData(4,:);
width = 1;
numBins = 50 / width;
[meanUtdLR, meanNN, lpUtdLR, upUtdLR, lpNN, upNN] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
for i = 1:numBins
    idx = width * (i - 1) < wL & wL <= width * i;
    meanUtdLR(2 * i - 1) = mean(iLossesUtdLR(idx));
    meanNN(2 * i - 1) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i - 1) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i - 1) = prctile(iLossesNN1(idx), 75);
    meanUtdLR(2 * i) = mean(iLossesUtdLR(idx));
    meanNN(2 * i) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i) = prctile(iLossesNN1(idx), 75);
    x(2 * i:2 * i + 1) = width * i;
end
x = x(1:end - 1);

%x = width:width:50;
figure
plot(x, meanNN)
hold on
grid on
plot(x, meanUtdLR, ':')
plot(x, lpNN)
plot(x, lpUtdLR, ':')
plot(x, upNN)
plot(x, upUtdLR, ':')
legend('NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR')
xlabel('Wedge length (m)')
ylabel('Mean absolute error (dBA)')
ylim([0 10])

%saveas(gcf, [savePath, filesep, 'WedgeLengthError'], 'epsc')
%saveas(gcf, [savePath, filesep, 'WedgeLengthError'], 'svg')

%% Function of minimum angle

%close all

[x, idx] = sort(inputData(3,:));
UtdLR = iLossesUtdLR(idx);
NN = iLossesNN2(idx);
y = [NN; UtdLR];

mA = inputData(3,:);
width = deg2rad(1);
numBins = deg2rad(90) / width;
[meanUtdLR, meanNN, lpUtdLR, upUtdLR, lpNN, upNN] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
num = zeros(1, numBins);
for i = 1:numBins
    idx = width * (i - 1) < mA & mA <= width * i;
    meanUtdLR(2 * i - 1) = mean(iLossesUtdLR(idx));
    meanNN(2 * i - 1) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i - 1) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i - 1) = prctile(iLossesNN1(idx), 75);
    meanUtdLR(2 * i) = mean(iLossesUtdLR(idx));
    meanNN(2 * i) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i) = prctile(iLossesNN1(idx), 75);
    x(2 * i:2 * i + 1) = width * i;
    num(i) = sum(idx);
end
x = x(1:end - 1);

%x = width:width:50;
figure
plot(x, meanNN)
hold on
grid on
plot(x, meanUtdLR, ':')
plot(x, lpNN)
plot(x, lpUtdLR, ':')
plot(x, upNN)
plot(x, upUtdLR, ':')
legend('NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR')
xlabel('Minimum Angle (rad)')
ylabel('Mean absolute error (dBA)')
ylim([0 10])

%saveas(gcf, [savePath, filesep, 'MinAngleError'], 'epsc')
%saveas(gcf, [savePath, filesep, 'MinAngleError'], 'svg')

%% Function of bending angle

%close all

[x, idx] = sort(inputData(2,:));
UtdLR = iLossesUtdLR(idx);
NN = iLossesNN1(idx);
y = [NN; UtdLR];

bA = inputData(2,:);
width = 0.1;
start = 3.1;
numBins = round(3.3 / width);
[meanUtdLR, meanNN, lpUtdLR, upUtdLR, lpNN, upNN] = deal(zeros(1, 2 * numBins));
x = start * ones(1, 2 * numBins + 1);
for i = 1:numBins
    idx = start + width * (i - 1) < bA & bA <= start + width * i;
    meanUtdLR(2 * i - 1) = mean(iLossesUtdLR(idx));
    meanNN(2 * i - 1) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i - 1) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i - 1) = prctile(iLossesNN1(idx), 75);
    meanUtdLR(2 * i) = mean(iLossesUtdLR(idx));
    meanNN(2 * i) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i) = prctile(iLossesNN1(idx), 75);
    x(2 * i:2 * i + 1) = start + width * i;
end
x = x(1:end - 1);

%x = width:width:50;
figure
plot(x, meanNN)
hold on
grid on
plot(x, meanUtdLR, ':')
plot(x, lpNN)
plot(x, lpUtdLR, ':')
plot(x, upNN)
plot(x, upUtdLR, ':')
legend('NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR')
xlabel('Bending Angle (rad)')
ylabel('Mean absolute error (dBA)')
ylim([0 10])

%saveas(gcf, [savePath, filesep, 'BendingAngleError'], 'epsc')
%saveas(gcf, [savePath, filesep, 'BendingAngleError'], 'svg')

%% Function of receiver angle

%close all

[x, idx] = sort(inputData(2,:) + inputData(3,:));
UtdLR = iLossesUtdLR(idx);
NN = iLossesNN1(idx);
y = [NN; UtdLR];

rA = inputData(2,:) + inputData(3,:);
width = 0.1;
start = 3.1;
numBins = round(3.3 / width);
[meanUtdLR, meanNN, lpUtdLR, upUtdLR, lpNN, upNN] = deal(zeros(1, 2 * numBins));
x = start * ones(1, 2 * numBins + 1);
for i = 1:numBins
    idx = start + width * (i - 1) < rA & rA <= start + width * i;
    meanUtdLR(2 * i - 1) = mean(iLossesUtdLR(idx));
    meanNN(2 * i - 1) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i - 1) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i - 1) = prctile(iLossesNN1(idx), 75);
    meanUtdLR(2 * i) = mean(iLossesUtdLR(idx));
    meanNN(2 * i) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i) = prctile(iLossesNN1(idx), 75);
    x(2 * i:2 * i + 1) = start + width * i;
end
x = x(1:end - 1);

%x = width:width:50;
figure
plot(x, meanNN)
hold on
grid on
plot(x, meanUtdLR, ':')
plot(x, lpNN)
plot(x, lpUtdLR, ':')
plot(x, upNN)
plot(x, upUtdLR, ':')
legend('NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR')
legend('NN-IIR (proposed)', 'UTD-LR')
xlabel('Receiver Angle (rad)')
ylabel('Mean absolute error (dBA)')
ylim([0 10])

%saveas(gcf, [savePath, filesep, 'RAngleError'], 'epsc')
%saveas(gcf, [savePath, filesep, 'RAngleError'], 'svg')

%% Function of path length

%close all

L = sqrt((inputData(5,:) + inputData(6,:)) .^ 2 + abs(inputData(7,:) - inputData(8,:)) .^ 2);

[x, idx] = sort(L);
UtdLR = iLossesUtdLR(idx);
NN = iLossesNN1(idx);
y = [NN; UtdLR];

pL = L;
width = 1;
numBins = 100 / width;
[meanUtdLR, meanNN, lpUtdLR, upUtdLR, lpNN, upNN] = deal(zeros(1, 2 * numBins));
x = ones(1, 2 * numBins + 1);
for i = 1:numBins
    idx = width * (i - 1) < pL & pL <= start + width * i;
    meanUtdLR(2 * i - 1) = mean(iLossesUtdLR(idx));
    meanNN(2 * i - 1) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i - 1) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i - 1) = prctile(iLossesNN1(idx), 75);
    meanUtdLR(2 * i) = mean(iLossesUtdLR(idx));
    meanNN(2 * i) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i) = prctile(iLossesNN1(idx), 75);
    x(2 * i:2 * i + 1) = width * i;
end
x = x(1:end - 1);

%x = width:width:50;
figure
plot(x, meanNN)
hold on
grid on
plot(x, meanUtdLR, ':')
plot(x, lpNN)
plot(x, lpUtdLR, ':')
plot(x, upNN)
plot(x, upUtdLR, ':')
ylim([0 5])
legend('NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR')
legend('NN-IIR (proposed)', 'UTD-LR')
xlabel('Path Length (m)')
ylabel('Mean absolute error (dBA)')
ylim([0 10])

%saveas(gcf, [savePath, filesep, 'PathLengthError'], 'epsc')
%saveas(gcf, [savePath, filesep, 'PathLengthError'], 'svg')

%% Function of wedge angle

%close all

[x, idx] = sort(inputData(1,:));
UtdLR = iLossesUtdLR(idx);
NN = iLossesNN1(idx);
y = [NN; UtdLR];

wA = inputData(1,:);
width = 0.1;
start = 3.1;
numBins = round(3.3 / width);
[meanUtdLR, meanNN, lpUtdLR, upUtdLR, lpNN, upNN] = deal(zeros(1, 2 * numBins));
x = start * ones(1, 2 * numBins + 1);
for i = 1:numBins
    idx = start + width * (i - 1) < wA & wA <= start + width * i;
    meanUtdLR(2 * i - 1) = mean(iLossesUtdLR(idx));
    meanNN(2 * i - 1) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i - 1) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i - 1) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i - 1) = prctile(iLossesNN1(idx), 75);
    meanUtdLR(2 * i) = mean(iLossesUtdLR(idx));
    meanNN(2 * i) = mean(iLossesNN1(idx));
    lpUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 25);
    upUtdLR(2 * i) = prctile(iLossesUtdLR(idx), 75);
    lpNN(2 * i) = prctile(iLossesNN1(idx), 25);
    upNN(2 * i) = prctile(iLossesNN1(idx), 75);
    x(2 * i:2 * i + 1) = start + width * i;
end
x = x(1:end - 1);

%x = width:width:50;
figure
plot(x, meanNN)
hold on
grid on
plot(x, meanUtdLR, ':')
plot(x, lpNN)
plot(x, lpUtdLR, ':')
plot(x, upNN)
plot(x, upUtdLR, ':')
legend('NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR', 'NN-IIR (proposed)', 'UTD-LR')
legend('NN-IIR (proposed)', 'UTD-LR')
xlabel('Wedge Angle (rad)')
ylabel('Mean absolute error (dBA)')
ylim([0 10])

%saveas(gcf, [savePath, filesep, 'WedgeAngleError'], 'epsc')
%saveas(gcf, [savePath, filesep, 'WedgeAngleError'], 'svg')

percentileUTDLRWedge = prctile(UtdLR, percentiles);
percentileNNWedge = prctile(NN, percentiles);

return

figure
plot(inputData(4,:), iLossesUtdLR, 'o')
hold on
plot(inputData(4,:), iLossesNN1, 'x')
legend('UTD-LR', 'NN-IIR (proposed)')

for i = 1:8
    figure
    plot(inputData(i,:), iLossesUtdLR, 'o')
    hold on
    plot(inputData(i,:), iLossesNN1, 'x')
end
