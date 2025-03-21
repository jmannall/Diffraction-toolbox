
close all
clear all
set(0, 'DefaultLineLineWidth', 1.5);
colorStore = colororder;

%% Constants

fs = 48e3;
nfft = 8192;
c = 344;
testSize = 20e3;
nBands = 8;

%% NN Parameters

loadDir = 'NNSaves';
netNames = {['Run6', filesep, 'IIR-7_32_0001.mat'], ['Run5', filesep, 'IIR-7_27_0001.mat'], ['Run6', filesep, 'IIR-5_45_0001.mat'], ['Run3', filesep, 'IIR-4_31_0001.mat'], ['Run6', filesep, 'IIR-4_20_0001.mat']};
netNames = {['Run6', filesep, 'IIR-7_32_0001.mat'], ['Run6', filesep, 'IIR-4_20_0001.mat'], 'iir-2057_0001-1-09-099-3-25.mat'};
netNames = {['Run6', filesep, 'IIR-7_32_0001.mat'], ['Run6', filesep, 'IIR-4_20_0001.mat']};

loadDir = 'NNSaves_FinalRun';
netNames = {['Run6', filesep, 'IIR-7_36_0001.mat'], ['Run6', filesep, 'IIR-4_20_0001.mat']};
% loadDir = 'NNSaves_UniformZw';
% netNames = {['Run2', filesep, 'IIR-7_27_0001.mat'], ['Run5', filesep, 'IIR-6_16_0001.mat']};
% loadDir = 'NNSaves_75Uniform_Zw';
% netNames = {['Run5', filesep, 'IIR-7_27_0001.mat'], ['Run5', filesep, 'IIR-4_20_0001.mat']};
% loadDir = 'NNSaves_Times2';
% netNames = {['Run2', filesep, 'IIR-4_52_0001.mat'], ['Run3', filesep, 'IIR-7_14_0001.mat']};

%loadDir = 'NNSaves_DC';
%netNames = {['Run2', filesep, 'IIR-6_40_0001.mat'], ['Run6', filesep, 'IIR-4_20_0001.mat']};

%% Geometry input data

disp('Load validation data')

controlparameters = struct('fs', 2 * fs, 'nfft', 2 * nfft, 'difforder', 1, 'c', c, 'saveFiles', 3, 'noDirect', false);
[inputData, targetData, ~, fvec, ~, ~, ~, ~, tfmag.Btm, tfmag.BtmI] = CreateBtmTrainingData(testSize, controlparameters, 'ValidationData');

numFreq = length(fvec);

controlparameters.fs = fs;
controlparameters.nfft = nfft;
[~, tfmagDefault, ~, fvec] = DefaultBTM(controlparameters);

[~, fc, fidx] = CreateFrequencyNBands(tfmagDefault, fvec, nBands);

wI = inputData(1,:);
bA = inputData(2,:);
mA = inputData(3,:);
wL = inputData(4,:);
rS = inputData(5,:);
rR = inputData(6,:);
zS = inputData(7,:);
zR = inputData(8,:);
X = dlarray(single(inputData), "CB");
clear inputData

wI = rad2deg(wI);
mA = rad2deg(mA);
bA = rad2deg(bA);
thetaS = mA;
thetaR = mA + bA;

[~, phii] = CalculateApex(rS, rR, zS, zR, wL, false);
dZ = abs(zR - zS);
pathLength = sqrt((rS + rR) .^ 2 + dZ .^ 2);
tfmag.Btm = mag2db((1 ./ pathLength) .* db2mag(tfmag.Btm));
tfmag.BtmI = mag2db((1 ./ pathLength) .* db2mag(tfmag.BtmI));

%% BTM

disp('BTM')

tfmagN.BtmI = mag2db((1 ./ pathLength) .* db2mag(targetData));
tfmagN.Btm = CreateFrequencyNBands(tfmag.Btm, fvec, nBands);

%% NN

disp('NN')

numFilters = 2;

numNets = length(netNames);
for i = 1:numNets
    disp(['Net ', num2str(i)])
    load([loadDir, filesep, netNames{i}], 'net');
    field = ['NN', num2str(i)];
    [tfmag.(field), tfmagN.(field)] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);
end

%% UTD LR

disp('UTD LR')

tf = zeros(testSize, 4);
if isfield(controlparameters, 'fvec')
    controlparameters = rmfield(controlparameters, 'fvec');
end
for i = 1:testSize
    tf(i,:) = SingleUTDWedgeInterpolated(thetaS(i), thetaR(i), rS(i), rR(i), wI(i), phii(i), controlparameters);
end

validPath = true(size(pathLength));
input = [1; zeros(5, 1)];
windowLength = length(input);

[~, ~, tfmag.UtdILR] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tf, c, fs, validPath);
[tfmagN.UtdILR, ~, ~] = CreateFrequencyNBands(tfmag.UtdILR, fvec, nBands);

controlparameters.fvec = fvec;
clear tf validPath input windowLength

%% UTD

disp('UTD')

fvecUtd = fvec(2:end);
controlparameters.fvec =  fvecUtd;
tfmag.UtdI = zeros(length(fvecUtd), testSize);
for i = 1:testSize
    tfmag.UtdI(:,i) = SingleUTDWedgeInterpolated(thetaS(i), thetaR(i), rS(i), rR(i), wI(i), phii(i), controlparameters);
end
[tfmagN.UtdI, ~, ~] = CreateFrequencyNBands(tfmag.UtdI, fvecUtd, nBands);
controlparameters.fvec = fvec;

%% IIR filter Kirsch and Ewert

disp('IIR')

tfmag.IIR = zeros(length(fvec), testSize);
for i = 1:testSize
    tfmag.IIR(:,i) = SingleIIRWedge(wL(i), wI(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), 4, controlparameters);
end
tfmagN.IIR = CreateFrequencyNBands(tfmag.IIR, fvec, nBands);

% disp('IIRHi')
% 
% n = 4;
% tfmag.IIRHi = zeros(length(fvec), testSize);
% for i = 1:testSize
%     tfmag.IIRHi(:,i) = SingleIIRWedge(wL(i), wI(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), n, controlparameters);
% end
% tfmagN.IIRHi = CreateFrequencyNBands(tfmag.IIRHi, fvec, nBands);
% 
% disp('IIRLo')
% 
% n = 2;
% tfmag.IIRLo = zeros(length(fvec), testSize);
% for i = 1:testSize
%     tfmag.IIRLo(:,i) = SingleIIRWedge(wL(i), wI(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), n, controlparameters);
% end
% tfmagN.IIRLo = CreateFrequencyNBands(tfmag.IIRLo, fvec, nBands);

%% Plots constants

savePath = 'figures/EURASIP';
fontSize = 18;
width = 300;
gap = 10;
numSets = 5;
fields = fieldnames(tfmag);
numFields = length(fields);
n = 3;
m = 4;

roundInput = double(round([rad2deg(extractdata(X(1:3,:))); extractdata(X(4:8,:))], 1));

%% Frequency Plots

disp('Frequency Plots')

close all

for j = 1:numSets
    figure('Position', [gap gap 6 * width 3 * width])
    t = tiledlayout(n ,m, 'TileSpacing', 'compact');
    title(t, 'Frequency responses for given input parameters')
    for i = 1:n * m
        idx = (j - 1) * n * m + i;
        nexttile
        for k = 1:numFields
            field = fields{k};
            if matches(field, 'UtdI')
                f = fvecUtd;
            else
                f = fvec;
            end
            semilogx(f, tfmag.(field)(:,idx))
            hold on
        end
        grid on
        l = legend(fields, 'Location','northwest');
        title(l, 'Diffraction Model')
        xlim([50 20e3])
        ylim([-70 0])
        ylabel('Magnitude (dB)')
        xlabel('Frequency (Hz)')
        title({['wL: ' num2str(roundInput(4, idx)), ' wI: ', num2str(roundInput(1, idx)), ' mA: ', num2str(roundInput(3, idx)), ' bA: ', num2str(roundInput(2, idx))] ...
            [' rS: ' num2str(roundInput(5, idx)), ' rR: ', num2str(roundInput(6, idx)), ' zS: ', num2str(roundInput(7, idx)), ' zR: ', num2str(roundInput(8, idx))]})
    end
    saveas(gcf, [savePath, filesep, 'AccuracyExample_', num2str(j), '.svg'], 'svg')
end

%% Losses

lossInf = CalculateLoss(tfmagN, tfmagN.BtmI, wL > 5);
lossN = CalculateLoss(tfmagN, tfmagN.BtmI);
lossNTrue = CalculateLoss(tfmagN, tfmagN.Btm);

% lossN.if.IIRHi = lossNTrue.if.IIRHi;
% lossN.i.IIRHi = mean(lossN.if.IIRHi, 1);
% lossN.f.IIRHi = mean(lossN.if.IIRHi, 2);
% lossN.mean.IIRHi = mean(lossN.if.IIRHi, 'all');
% 
% lossN.if.IIRLo = lossNTrue.if.IIRLo;
% lossN.i.IIRLo = mean(lossN.if.IIRLo, 1);
% lossN.f.IIRLo = mean(lossN.if.IIRLo, 2);
% lossN.mean.IIRLo = mean(lossN.if.IIRLo, 'all');

%% Frequency dependent error

%close all

figure
for k = 3:numFields
    field = fields{k};
    if matches(field, 'UtdI')
        f = fvecUtd;
    else
        f = fvec;
    end
    semilogx(fc, lossInf.f.(field))
    hold on
end
grid on
l = legend(fields{3:numFields}, 'Location','northwest');
title(l, 'Diffraction Model')
xlim([20 20e3])
ylim([0 12])

figure
for k = 3:numFields
    field = fields{k};
    if matches(field, 'UtdI')
        f = fvecUtd;
    else
        f = fvec;
    end
    semilogx(fc, lossN.f.(field))
    hold on
end
grid on
l = legend(fields{3:numFields}, 'Location','northwest');
title(l, 'Diffraction Model')
xlim([20 20e3])
ylim([0 12])

%% Histograms

[~,edges] = histcounts(lossN.i.UtdILR);

limit = [0 12e3];
for i = 2:numFields
    field = fields{i};
    figure
    histogram(lossN.i.(field), edges)
    title([field, ' losses'])
    xlabel('L1 norm loss')
    ylabel('Number of cases')
    ylim(limit)
    grid on
end

%saveas(gcf, [savePath, filesep, 'AccuracyHistograms.svg'], 'svg')

%% Shaded plots

percentiles = 10:10:50;
wLMin = 5;

idx = find(wL > 4.5 & wL < 5.5 & phii > 85 & rS < 3);
num = length(idx);

for i = 1:num

    figure
    semilogx(fvec, tfmag.BtmI(:,idx(i)))
    hold on
    semilogx(fvec, tfmag.NN1(:,idx(i)))
    semilogx(fvec, tfmag.UtdILR(:,idx(i)))
    ylim([-60 10])
    xlim([20 20e3])
    grid on
    legend('BtmI', 'NN', 'UtdILR')
end
close all
for i = 2:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field), wL, 1, percentiles, ['wL: ', field], i)
end

%%

[zA, phii] = CalculateApex(rS, rR, zS, zR, wL, true);
zA = zA(:,3);
%zA = zA ./ wL';

distance = min(zA, wL' - zA)';

idx = wL > 4 & wL < 6 & distance < 1;
fields = {'NN1', 'UtdILR'};
numFields = 2;
percentiles = 10:10:50;

close all
for i = 1:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(idx), phii(idx), 1, percentiles, ['phii: ', field], i)
end

idx = wL > 4 & wL < 6;
for i = 1:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(idx), rS(idx), 1, percentiles, ['rS: ', field], i)
end

for i = 1:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(idx), rR(idx), 1, percentiles, ['rR: ', field], i)
end

close all
for i = 2:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(wL > wLMin), mA(wL > wLMin), 1, percentiles, ['mA: ', field], i)
end
%%

close all
for i = 2:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(wL > wLMin), bA(wL > wLMin), 2, percentiles, ['bA: ', field], i)
end
%%

close all
for i = 2:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(wL > wLMin), wI(wL > wLMin), 2, percentiles, ['wI: ', field], i)
end

close all
for i = 2:numFields
    field = fields{i};
    figure
    PlotVarDependentLoss(lossN.i.(field)(wL > wLMin), pathLength(wL > wLMin), 1, percentiles, ['L: ', field], i)
end

%% Percentiles

close all

percentiles = 0:0.1:100;

percentile = CalculatePercentiles(lossN.i, percentiles);
percentileTrue = CalculatePercentiles(lossNTrue.i, percentiles);

color = colorStore([1, 1, 2, 5], :);

figure
colororder(color)
% for i = 2:numFields
%     idx = fields{i};
%     plot(percentiles, percentile.(idx))
%     hold on
% end
plot(percentiles, percentile.NN1)
hold on
plot(percentiles, percentile.NN2, 'Color', [color(1,:), 0.6])
plot(percentiles, percentile.UtdILR, '--')
%plot(percentiles, percentile.UtdI, '-.')
%plot(percentiles, percentile.IIRHi, ':')
plot(percentiles, percentileTrue.IIR, ':')
grid on
l = legend('NN (best)', 'NN (small)', 'UTD-LR', 'UTD-IIR', 'Location','northwest');
title(l, 'Diffraction Model')
xlabel('Percentile')
ylabel('Mean absolute error \Psi (dB)')
ylim([0 30])

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.2,0.26,0.3,0.35])
box(ax,'on')
% for i = 2:numFields
%     idx = fields{i};
%     plot(percentiles, percentile.(idx), 'parent',ax)
%     hold on
% end
plot(percentiles, percentile.NN1)
hold on
plot(percentiles, percentile.NN2, 'Color', [color(1,:), 0.6])
plot(percentiles, percentile.UtdILR, '--')
%plot(percentiles, percentile.UtdI, '-.')
%plot(percentiles, percentile.IIRHi, ':')
plot(percentiles, percentileTrue.IIR, ':')
grid on
set(ax,'xlim',[0,40],'ylim',[0,1])

saveas(gcf, [savePath, filesep, 'PercentileLoss'], 'epsc')
saveas(gcf, [savePath, filesep, 'PercentileLoss'], 'svg')

figure
for i = 3:numFields
    idx = fields{i};
    plot(percentiles, percentile.(idx))
    hold on
end
grid on
l = legend(fields{3:end}, 'Location','northwest');
title(l, 'Diffraction Model')
xlabel('Percentile')
ylabel('Mean absolute error \Psi (dB)')
ylim([0 40])

%saveas(gcf, [savePath, filesep, 'PercentileLoss_Cropped.svg'], 'svg')

%% Score

%close all

% score = sum(iLossesNN2 < iLossesUtdLR) / testSize * 100;
% comparison = iLossesUtdLR - iLossesNN2;
% 
% figure
% bar(comparison)
% grid on
% title('Improvement in NN-IIR over UTD-LR')
% ylabel('Mean absolute error improvement (dBA)')
% xlabel('Number of cases')
% ylim([-10 35])
% xticks(0:2e3:2e4)
% xticklabels(0:2e3:2e4)
% 
% %saveas(gcf, [savePath, filesep, 'Comparison.svg'], 'svg')
% 
% percentilesComparison = prctile(comparison, percentiles);
% 
% figure
% plot(percentiles, percentilesComparison)
% grid on
% xlabel('Mean absolute error improvement (dBA)')
% ylabel('Percentile')
% 
% figure
% histogram(comparison)
% grid on
% %title('Improvement in NN-IIR over UTD-LR')
% xlabel('Mean absolute error improvement (dBA)')
% ylabel('Number of cases')
% xlim([-10 35])

%saveas(gcf, [savePath, filesep, 'ComparisonHistogram'], 'epsc')
%saveas(gcf, [savePath, filesep, 'ComparisonHistogram'], 'svg')

%% Function of wedgeLength

close all

width = 1;
numBins = 50 / width;
[meanUtdLR, meanNN1, meanNN2, meanUTDIIR] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
for i = 1:numBins
    idx = width * (i - 1) < wL & wL <= width * i;
    meanUtdLR(2 * i - 1:2 * i) = mean(lossN.i.UtdILR(idx));
    meanNN1(2 * i - 1:2 * i) = mean(lossN.i.NN1(idx));
    meanNN2(2 * i - 1:2 * i) = mean(lossN.i.NN2(idx));
    meanUTDIIR(2 * i - 1:2 * i) = mean(lossNTrue.i.IIR(idx));
    x(2 * i:2 * i + 1) = width * i;
end
x = x(1:end - 1);

color = colorStore([1, 1, 2, 5], :);

%x = width:width:50;
figure
colororder(color)
plot(x, meanNN1)
hold on
grid on
plot(x, meanNN2, 'Color', [color(1,:), 0.6])
plot(x, meanUtdLR, '--')
plot(x, meanUTDIIR, ':')
%plot(x, lpNN)
%plot(x, lpUtdLR, ':')
%plot(x, upNN)
%plot(x, upUtdLR, ':')
%legend('NN-IIR (best)', 'NN-IIR (small)', 'UTD-LR', 'UTD-IIR')
xlabel('Wedge length (m)')
ylabel('Mean absolute error \Psi (dB)')
ylim([0 16])
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, 'WedgeLengthError'], 'epsc')
saveas(gcf, [savePath, filesep, 'WedgeLengthError'], 'svg')

%% Function of phii

close all

width = 2;
numBins = 90 / width;
[meanUtdLR, meanNN1, meanNN2, meanUTDIIR] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
for i = 1:numBins
    idx = width * (i - 1) < phii & phii <= width * i & wL > 4;
    meanUtdLR(2 * i - 1:2 * i) = mean(lossN.i.UtdILR(idx));
    meanNN1(2 * i - 1:2 * i) = mean(lossN.i.NN1(idx));
    meanNN2(2 * i - 1:2 * i) = mean(lossN.i.NN2(idx));
    meanUTDIIR(2 * i - 1:2 * i) = mean(lossNTrue.i.IIR(idx)); 
    x(2 * i:2 * i + 1) = width * i;
end
x = x(1:end - 1);

color = colorStore([1, 1, 2, 5], :);

%x = width:width:50;
figure
colororder(color)
plot(x, meanNN1)
hold on
grid on
plot(x, meanNN2, 'Color', [color(1,:), 0.6])
plot(x, meanUtdLR, '--')
plot(x, meanUTDIIR, ':')
%plot(x, lpNN)
%plot(x, lpUtdLR, ':')
%plot(x, upNN)
%plot(x, upUtdLR, ':')
legend('NN-IIR (best)', 'NN-IIR (small)', 'UTD-LR', 'UDFA')
xlabel('Phi (degrees)')
ylabel('Mean absolute error \Psi (dB)')
ylim([0 7])
xlim([0 90])
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, 'PhiiError'], 'epsc')
saveas(gcf, [savePath, filesep, 'PhiiError'], 'svg')

%% function of distance from corner

close all

[zA, phii] = CalculateApex(rS, rR, zS, zR, wL, true);
zA = zA(:,3);
%zA = zA ./ wL';

distance = min(zA, wL' - zA)';

width = 1;
numBins = 25 / width;
[meanUtdLR, meanNN1, meanNN2, meanUTDIIR] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
for i = 1:numBins
    idx = width * (i - 1) < distance & distance <= width * i & wL > 4;
    meanUtdLR(2 * i - 1:2 * i) = mean(lossN.i.UtdILR(idx));
    meanNN1(2 * i - 1:2 * i) = mean(lossN.i.NN1(idx));
    meanNN2(2 * i - 1:2 * i) = mean(lossN.i.NN2(idx));
    meanUTDIIR(2 * i - 1:2 * i) = mean(lossNTrue.i.IIR(idx)); 
    x(2 * i:2 * i + 1) = width * i;
end
x = x(1:end - 1);

color = colorStore([1, 1, 2, 5], :);

%x = width:width:50;
figure
colororder(color)
plot(x, meanNN1)
hold on
grid on
plot(x, meanNN2, 'Color', [color(1,:), 0.6])
plot(x, meanUtdLR, '--')
plot(x, meanUTDIIR, ':')
%plot(x, lpNN)
%plot(x, lpUtdLR, ':')
%plot(x, upNN)
%plot(x, upUtdLR, ':')
%legend('NN-IIR (best)', 'NN-IIR (small)', 'UTD-LR', 'UTD-IIR')
xlabel('Minimum distance of z_a from z_0 or z_w (m)')
ylabel('Mean absolute error \Psi (dB)')
ylim([0 4])
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, 'ZaError'], 'epsc')
saveas(gcf, [savePath, filesep, 'ZaError'], 'svg')

return

%% function of zA

close all

[zA, phii] = CalculateApex(rS, rR, zS, zR, wL, true);
zA = zA(:,3);
%zA = zA ./ wL';

zA = min(zA, wL' - zA);

percentiles = 10:10:50;

idx = wL > 4;
figure
PlotVarDependentLoss(lossN.i.NN1(idx), zA(idx), 0.5, percentiles, 'zA: NN', 1)

figure
PlotVarDependentLoss(lossN.i.UtdILR(idx), zA(idx), 0.5, percentiles, 'zA: UTD', 1)

idx = wL > 4 & phii > 45;

figure
PlotVarDependentLoss(lossN.i.UtdILR(idx), zA(idx), 0.5, percentiles, 'zA (phii > 45): UTD', 1)
idx = wL > 4 & phii > 80;

figure
PlotVarDependentLoss(lossN.i.UtdILR(idx), zA(idx), 0.5, percentiles, 'zA (phii > 80): UTD', 1)


idx = wL > 4;

figure
PlotVarDependentLoss(lossN.i.NN1(idx), phii(idx), 1, percentiles, 'phii: NN', 1)

figure
PlotVarDependentLoss(lossN.i.UtdILR(idx), phii(idx), 1, percentiles, 'phii: UTD', 1)


width = 0.1;
numBins = 1 / width;
[meanUtdLR, meanNN1, meanNN2, meanUTDIIR] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
for i = 1:numBins
    idx = width * (i - 1) < zA & zA <= width * i;
    meanUtdLR(2 * i - 1:2 * i) = mean(lossN.i.UtdILR(idx));
    meanNN1(2 * i - 1:2 * i) = mean(lossN.i.NN1(idx));
    meanNN2(2 * i - 1:2 * i) = mean(lossN.i.NN2(idx));
    meanUTDIIR(2 * i - 1:2 * i) = mean(lossNTrue.i.IIR(idx));
    x(2 * i:2 * i + 1) = width * i;
end
x = x(1:end - 1);

color = colorStore([1, 1, 2, 5], :);

%x = width:width:50;
figure
colororder(color)
plot(x, meanNN1)
hold on
grid on
plot(x, meanNN2, 'Color', [color(1,:), 0.6])
plot(x, meanUtdLR, '--')
plot(x, meanUTDIIR, ':')
%plot(x, lpNN)
%plot(x, lpUtdLR, ':')
%plot(x, upNN)
%plot(x, upUtdLR, ':')
legend('NN-IIR (best)', 'NN-IIR (small)', 'UTD-LR', 'UTD-IIR')
xlabel('zA')
ylabel('Mean absolute error (dB)')
ylim([0 10])

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
