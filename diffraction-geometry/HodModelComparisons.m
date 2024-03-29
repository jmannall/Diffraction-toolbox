%close all
clear all

set(0, 'DefaultLineLineWidth', 1.5);
colorStore = colororder('default');

%% Data
fs = 48e3;
nfft = 8192;
c = 344;
numEdges = 2;

height = 5;
%height = 2.4;
[zS, zR] = deal(height / 2);
%[zS, zR] = deal(0.5);

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges, 'c', c, 'saveFiles', 3, 'noDirect', true, 'interpolated', false);

%% Generate geometry and controls

numPaths = 500;

index = DataHash({numPaths, fs, nfft, numEdges, height});
index = '2180cfff4ddad1de505214480a4535ce';
[loadPath, savePath] = deal(['geometry/NthOrderPaths_01mTo3m_', num2str(index), '.mat']);
restart = true;
generate = false;
plotFigures = false;
createPlot = false;
if restart
    i = 1;
end
n = i;

dtemplate = struct('rS', [], 'rR', [], 'W', [], 'L', [], 'thetaS', [], 'thetaR', [], 'wedgeIndex', []);
if generate && n == 1
    data = repmat(dtemplate, numPaths, 1);
elseif ~generate
    load(loadPath, 'data')
end

[~, tfmagDefault, ~, fvec] = DefaultBTM(controlparameters);
nBand = 8;

[~, fc, fidx] = CreateFrequencyNBands(tfmagDefault, fvec, nBand); 

%% NN parameters

loadDir = 'NNSaves';
netName = ['Run6', filesep, 'IIR-7_32_0001.mat'];

loadDir = 'NNSaves_FinalRun';
netName = ['Run6', filesep, 'IIR-7_36_0001.mat'];

% loadDir = 'NNSaves_DC';
% netName = ['Run2', filesep, 'IIR-6_40_0001.mat'];
% 
% loadDir = 'NNSaves_UniformZw';
% netName = ['Run2', filesep, 'IIR-7_27_0001.mat'];
% 
% loadDir = 'NNSaves_75Uniform_Zw';
% netName = ['Run5', filesep, 'IIR-7_27_0001.mat'];

% loadDir = 'NNSaves_Run2';
% netName = ['Run4', filesep, 'IIR-6_30_0001.mat'];

numFilters = 2;
pathLength = ones(1, numEdges);
load([loadDir, filesep, netName], 'net');

%% Geometry parameters

wL = [height, height];
epsilon = 1e-2;

%% UTD-LR parameters

input = [1; zeros(5, 1)];
windowLength = length(input);
validPath = true(1, numEdges);
phii = 90;

disp('Start')
for i = n:numPaths
    if generate
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, data(i)] = GenerateNthOrderPath(numEdges, height);

        wI = data(i).wedgeIndex;
        thetaS = data(i).thetaS;
        thetaR = data(i).thetaR;
        rS = data(i).rS;
        rR = data(i).rR; 
        W = data(i).W;
        data(i).L = rS + sum(W) + rR;
        mA = deg2rad(epsilon) * ones(1, numEdges);
        bA = deg2rad([wI(1) - thetaS, wI(2:numEdges - 1), thetaR]);
    else
        wI = data(i).wedgeIndex;
        thetaS = data(i).thetaS;
        thetaR = data(i).thetaR;
        rS = data(i).rS;
        rR = data(i).rR;
        W = data(i).W;
        data(i).L = rS + sum(W) + rR;
        mA = deg2rad(epsilon) * ones(1, numEdges);
        bA = deg2rad([wI(1) - thetaS, wI(2:numEdges - 1), thetaR]);

        [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid, vReceiver] = CreateNthOrderPathData(wI, thetaS, thetaR, rS, rR, W, height);
    end

    % Expand to include z variation. Requires calculating all the apex points.
    
    %% 2nd order BTM

    controlparameters.fs = 2 * fs;
    controlparameters.nfft = 2 * nfft;
    [~, tfmagStore, ~, fvecBtm, ~] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);
    tfmag.Btm = tfmagStore.diff2;
    tfmagN.Btm(:,i) = CreateNBandMagnitude(tfmag.Btm, fidx);

    %% BTM Daisy Chains
    
    % BTM with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [tfmag.BtmA, ~, tf.BtmA] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    tfmagN.BtmA(:,i) = CreateNBandMagnitude(tfmag.BtmA(:,end), fidx);

    % BTM with virtual sources and receivers set at rS + W and W + rR from the
    % edge and normalised
    [tfmag.BtmE, ~, tf.BtmE] = SingleBTMExtDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    tfmagN.BtmE(:,i) = CreateNBandMagnitude(tfmag.BtmE(:,end), fidx);

    %% BTM Daisy Chains interpolated

    [tfmag.BtmIE, tfmag.BtmIA, tf.BtmIE, tf.BtmIA] = deal(zeros(nfft, numEdges + 1));

    disp('BTM Ext')

    [tfmag.BtmIE(:,1), ~, tf.BtmIE(:,1)] = SingleWedgeInterpolated(wL(1), wI(1), epsilon, wI(1) - thetaS, rS, W + rR, zS, zR, controlparameters, false);
    [tfmag.BtmIE(:,2), ~, tf.BtmIE(:,2)] = SingleWedgeInterpolated(wL(2), wI(2), epsilon, thetaR, rS + W, rR, zS, zR, controlparameters, false);
    
    tf.BtmIE(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(tf.BtmIE(:,1:numEdges), 2);
    tfmag.BtmIE(:,end) = mag2db(abs(tf.BtmIE(:,end)));

    tfmagN.BtmIE(:,i) = CreateNBandMagnitude(tfmag.BtmIE(:,end), fidx);

    disp('BTM Apex')

    [tfmag.BtmIA(:,1), ~, tf.BtmIA(:,1)] = SingleWedgeInterpolated(wL(1), wI(1), epsilon, wI(1) - thetaS, rS, W, zS, zR, controlparameters, false);
    [tfmag.BtmIA(:,2), ~, tf.BtmIA(:,2)] = SingleWedgeInterpolated(wL(2), wI(2), epsilon, thetaR, W, rR, zS, zR, controlparameters, false);

    tf.BtmIA(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(tf.BtmIA(:,1:numEdges), 2);
    tfmag.BtmIA(:,end) = mag2db(abs(tf.BtmIA(:,end)));

    tfmagN.BtmIA(:,i) = CreateNBandMagnitude(tfmag.BtmIA(:,end), fidx);

    %% Neural networks

    controlparameters.fs = fs;
    controlparameters.nfft = nfft;
    [tfmag.NNE, tfmag.NNA, tf.NNE, tf.NNA] = deal(zeros(nfft / 2, numEdges + 1));

    disp('NN Ext')

    r1 = min([rS, rS + W], [W + rR, rR]);
    r2 = max([rS, rS + W], [W + rR, rR]);
    X = [deg2rad(wI); bA; mA; wL; r1; r2; [zS, zS]; [zR, zR]];
    %X = [deg2rad(wI); bA; mA; wL / 2; r1; r2; [0.01, 0.01]; [0.01, 0.01]];
    X = dlarray(single(X), "CB");

    [tfmag.NNE(:,1:numEdges), ~, tf.NNE(:,1:numEdges)] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);

    tf.NNE(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(tf.NNE(:,1:numEdges), 2);
    tfmag.NNE(:,end) = mag2db(abs(tf.NNE(:,end)));

    tfmagN.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,end), fidx);

    disp('NN Apex')

    r1 = min([rS, W], [W, rR]);
    r2 = max([rS, W], [W, rR]);
    X = [deg2rad(wI); bA; mA; wL; r1; r2; [zS, zS]; [zR, zR]];
    X = dlarray(single(X), "CB");

    [tfmag.NNA(:,1:numEdges), ~, tf.NNA(:,1:numEdges)] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);

    tf.NNA(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(tf.NNA(:,1:numEdges), 2);
    tfmag.NNA(:,end) = mag2db(abs(tf.NNA(:,end)));

    tfmagN.NNA(:,i) = CreateNBandMagnitude(tfmag.NNA(:,end), fidx);

    %% UTD
    
    controlparameters.interpolated = false;
    
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    tfmagStore = SingleUTDApexDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r

    [~, ~, tfmag.UtdLRA] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tfmagStore', c, fs, validPath);
    [tfmagN.UtdLRA(:,i), ~, ~] = CreateFrequencyNBands(tfmag.UtdLRA(:,end), fvec, nBand);

    tfmagStore = SingleUTDExtDaisyChain(data(i), phii, controlparameters);

    [~, ~, tfmag.UtdLRE] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tfmagStore', c, fs, validPath);
    [tfmagN.UtdLRE(:,i), ~, ~] = CreateFrequencyNBands(tfmag.UtdLRE(:,end), fvec, nBand);
    
    %% UTD Interpolated

    controlparameters.interpolated = true;
        
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    tfmagStore = SingleUTDApexDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r

    [~, ~, tfmag.UtdILRA] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tfmagStore', c, fs, validPath);
    [tfmagN.UtdILRA(:,i), ~, ~] = CreateFrequencyNBands(tfmag.UtdILRA(:,end), fvec, nBand);
    
    tfmagStore = SingleUTDExtDaisyChain(data(i), phii, controlparameters);

    [~, ~, tfmag.UtdILRE] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tfmagStore', c, fs, validPath);
    [tfmagN.UtdILRE(:,i), ~, ~] = CreateFrequencyNBands(tfmag.UtdILRE(:,end), fvec, nBand);
    
    %% Figures
    if plotFigures
        close all

        figure('Position',[50 100 900 800])
        plot3(source(:,1), source(:,2), source(:,3), 'o')
        hold on
        plot3(receiver(:,1), receiver(:,2), receiver(:,3), 'o')
        plot3(Q(:,1), Q(:,2), Q(:,3))
        plot3(apex(:,1), apex(:,2), apex(:,3), 'o')
        title('Scene')
        legend('Source', 'Receiver', 'Planes', 'Apex')
        view([0 90])
        xlim([-4 2])
        ylim([-4 2])

        color = colorStore(1:4,:);
        color = colorStore([4, 1, 7, 2],:);
        figure('Position',[50 100 1820 800])
        %tiledlayout(1,2);

        %nexttile
        semilogx(fvecBtm, tfmag.Btm, 'LineWidth', 2)
        hold on
        %semilogx(fvecBtm, tfmag.BtmE(:,end))
        %semilogx(fvec, tfmag.UtdLRE(:,end))
        semilogx(fvec, tfmag.NNE(:,end), '--')
        semilogx(fvecBtm, tfmag.BtmIE(:,end), '--')
        semilogx(fvec, tfmag.UtdILRE(:,end), '--')
        %semilogx(fvecBtm, tfmag.BtmIA(:,end), '-.')
        %semilogx(fvec, tfmag.UtdLRA(:,end), '-.')
        %semilogx(fvec, tfmag.NNA(:,end), ':')
        %semilogx(fvecBtm, tfmag.BtmIA(:,end), ':')
        %semilogx(fvec, tfmag.UtdILRA(:,end), ':')
        grid on
        colororder(gca, color)
        title('Frequency responses')
        xlabel('Frequency')
        ylabel('Magnitude')
        %legend('True BTM', 'BTM Ext', 'UTD Ext', 'NN Ext', 'BTMI Ext', 'UTDI Ext', 'BTM Apex', 'UTD Apex', 'NN Apex', 'BTM ApexI', 'UTD ApexI')
        legend('True BTM', 'NN Ext', 'BTMI Ext', 'UTDI Ext')
        xlim([20 20000])
        ylim([-70 0])

%             color = colorStore(1:5,:);
%             nexttile
%             for j = 1:2
%                 semilogx(fvec, tfmagNNExt(:,j), 'LineWidth', 2)
%                 hold on
%                 semilogx(fvecBtm, btmTfmagExt(:,j))
%                 semilogx(fvecBtm, btmTfmagExtI(:,j), '--')
%                 semilogx(fvec, utdTfmagExtLR(:,j), '-.')
%                 semilogx(fvec, utdTfmagExtILR(:,j), ':')
%             end
%             grid on
%             colororder(gca, color)
%             title('Frequency responses')
%             xlabel('Frequency')
%             ylabel('Magnitude')
%             legend('NN', 'BTM Ext', 'BTM ExtI', 'UTD Ext', 'UTD ExtI')
%             xlim([20 20000])
%             ylim([-30 10])

        color = colorStore([4, 4, 1, 1, 2, 2],:);
        figure
        %semilogx(fvecBtm, tfmag.BtmE(:,end))
        %semilogx(fvec, tfmag.UtdLRE(:,end))
        semilogx(fvecBtm, tfmag.BtmIE(:,1:2), 'LineWidth', 2)
        hold on
        semilogx(fvec, tfmag.NNE(:,1:2), '--')
        semilogx(fvec, tfmag.UtdILRE(:,1:2), '--')
        %semilogx(fvecBtm, tfmag.BtmIA(:,end), '-.')
        %semilogx(fvec, tfmag.UtdLRA(:,end), '-.')
        %semilogx(fvec, tfmag.NNA(:,end), ':')
        %semilogx(fvecBtm, tfmag.BtmIA(:,end), ':')
        %semilogx(fvec, tfmag.UtdILRA(:,end), ':')
        grid on
        colororder(color)
        title('Frequency responses')
        xlabel('Frequency')
        ylabel('Magnitude')
        %legend('True BTM', 'BTM Ext', 'UTD Ext', 'NN Ext', 'BTMI Ext', 'UTDI Ext', 'BTM Apex', 'UTD Apex', 'NN Apex', 'BTM ApexI', 'UTD ApexI')
        legend('BTMI Ext', '', 'NN Ext', '', 'UTDI Ext', '')
        xlim([20 20000])
        ylim([-60 10])
    end
    tfmagN1.BtmIE(:,i) = CreateNBandMagnitude(tfmag.BtmIE(:,1), fidx);
    tfmagN1.UtdILRE(:,i) = CreateFrequencyNBands(tfmag.UtdILRE(:,1), fvec, nBand);
    tfmagN1.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,1), fidx);

    tfmagN2.BtmIE(:,i) = CreateNBandMagnitude(tfmag.BtmIE(:,2), fidx);
    tfmagN2.UtdILRE(:,i) = CreateFrequencyNBands(tfmag.UtdILRE(:,2), fvec, nBand);
    tfmagN2.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,2), fidx);
end

%%
lossN1 = CalculateLoss(tfmagN1, tfmagN1.BtmIE);
lossN2 = CalculateLoss(tfmagN2, tfmagN2.BtmIE);

loss1.NNE = (lossN1.mean.NNE + lossN2.mean.NNE) / 2;
loss1.UtdILRE = (lossN1.mean.UtdILRE + lossN2.mean.UtdILRE) / 2;

loss = CalculateLoss(tfmagN, tfmagN.Btm);

% for i = 1:numPaths
%     disp(['NN Loss: ', num2str(lossN1.i.NNE), ' + ', num2str(lossN2.i.NNE), ' = ', num2str(loss.i.NNE(i))])
%     disp(['UTD Loss: ', num2str(lossN1.i.UtdILRE), ' + ', num2str(lossN2.i.UtdILRE), ' = ', num2str(loss.i.UtdILRE(i))])
% end
disp('Complete')
    
%% Data

%close all

if isfield(loss, 'w')
    loss = rmfield(loss, 'w');
end

L = [data.L];
W = [data.W];
W = W(1:numPaths);
fields = fieldnames(loss.i);
numFields = length(fields);

width = 0.1;
numBins = 3 / width;
x = (0:width:3 + width) - width / 2;
for i = 2:numBins+1
    idx = width * (i - 2) < W & W <= width * (i - 1);
    for j = 1:numFields
        field = fields{j};
        loss.w.(field)(i) = mean(loss.i.(field)(idx));
    end
end
for j = 1:numFields
    field = fields{j};
    idx = isnan(loss.w.(field));
    loss.w.(field)(idx) = 0;
    loss.w.(field)(numBins + 2) = 0;
end

%% Plot
saveDir = 'figures';
color = colorStore([4,4,1,2,4,4,1,2],:);
f = figure;
plot(x, loss.w.BtmE, 'Color', [color(1,:), 0.6])
hold on
grid on
plot(x, loss.w.BtmIE)
plot(x, loss.w.NNE)
plot(x, loss.w.UtdILRE)
plot(x, loss.w.BtmA, '--', 'Color', [color(1,:), 0.6])
plot(x, loss.w.BtmIA, '--')
plot(x, loss.w.NNA, '--')
plot(x, loss.w.UtdILRA, '--')

colororder(color)
xlim([0.1 3])
ylim([0 6])
legend('BTM Extension', 'BTM-I Extension', 'NN-IIR (best) Extension', 'UTD-LR Extension', 'BTM Apex', 'BTM-I Apex', 'NN-IIR (best) Apex', 'UTD-LR Apex')
xlabel('W_1 (m)')
ylabel('Mean absolute error \Psi (dB)')
%fontsize(f, 16, "points")

saveas(gcf, [saveDir filesep 'HODComparison'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparison'], 'svg')

%%

if isfield(loss, 'tR')
    loss = rmfield(loss, 'tR');
end

for i = 1:numPaths
thetaR(i) = data(i).wedgeIndex(1);
end
thetaR = [data.L];

width = 180 / 20;
numBins = 180 / width;
x = (180:width:360 + width) - width / 2;
width = 0.2;
numBins = 5 / width;
x = (0:width:5 + width) - width / 2;

for i = 2:numBins+1
    idx = x(i - 1) < thetaR & thetaR <= x(i);
    for j = 1:numFields
        field = fields{j};
        loss.tR.(field)(i) = mean(loss.i.(field)(idx));
    end
end
for j = 1:numFields
    field = fields{j};
    idx = isnan(loss.w.(field));
    loss.tR.(field)(idx) = 0;
    loss.tR.(field)(numBins + 2) = 0;
end

saveDir = 'figures';
color = colorStore([4,4,1,2,5,5,6,3],:);

figure
plot(x, loss.tR.BtmE)
hold on
grid on
plot(x, loss.tR.BtmIE, 'Color', [color(1,:), 0.6])
plot(x, loss.tR.NNE)
plot(x, loss.tR.UtdILRE)
plot(x, loss.tR.BtmA, '--')
plot(x, loss.tR.BtmIA, '--', 'Color', [color(5,:), 0.6])
plot(x, loss.tR.NNA, '--')
plot(x, loss.tR.UtdILRA, '--')
colororder(color)
%xlim([180 360])
%ylim([0 6])
legend('BTM Extension', 'BTM-I Extension', 'NN-IIR Extension', 'UTD-LR Extension', 'BTM Apex', 'BTM-I Apex', 'NN-IIR Apex', 'UTD-LR Apex')
xlabel('W_1 (m)')
ylabel('Mean absolute error (dB)')

%%

color = colorStore([4,4,1,2,5,5,6,3],:);

figure
semilogx(fc, loss.f.BtmE)
hold on
grid on
semilogx(fc, loss.f.BtmIE, 'Color', [color(1,:), 0.6])
semilogx(fc, loss.f.NNE)
semilogx(fc, loss.f.UtdILRE)
semilogx(fc, loss.f.BtmA, '--')
semilogx(fc, loss.f.BtmIA, '--', 'Color', [color(5,:), 0.6])
semilogx(fc, loss.f.NNA, '--')
semilogx(fc, loss.f.UtdILRA, '--')
colororder(color)
legend('BTM Extension', 'BTM-I Extension', 'NN-IIR Extension', 'UTD-LR Extension', 'BTM Apex', 'BTM-I Apex', 'NN-IIR Apex', 'UTD-LR Apex')
xlim([62.5 20e3])
grid on

color = colorStore([1,2],:);

figure
semilogx(fc, lossN1.f.NNE)
hold on
grid on
semilogx(fc, lossN1.f.UtdILRE)
colororder(color)
legend('NN-IIR Extension', 'UTD-LR Extension')
xlim([62.5 20e3])
grid on

figure
semilogx(fc, lossN2.f.NNE)
hold on
grid on
semilogx(fc, lossN2.f.UtdILRE)
colororder(color)
legend('NN-IIR Extension', 'UTD-LR Extension')
xlim([62.5 20e3])
grid on

return

close all
saveDir = 'figures';

for k = 1:length(store)
    load(['hod workspace_Test_' num2str(k), '.mat'])
    
    color = colorStore([4,1,2,5,6,3],:);

    figure
    plot(x, meanBE)
    hold on
    grid on
    plot(x, meanNE)
    plot(x, meanUE)
    plot(x, meanBA, '--')
    plot(x, meanNA, '--')
    plot(x, meanUA, '--')
    %plot(x, meanUEI)
    colororder(color)
    xlim([0.1 3])
    ylim([0 5])
    title(num2str(store(k)))
    legend('BTM Extension', 'NN-IIR Extension', 'UTD-LR Extension', 'BTM Apex', 'NN-IIR Apex', 'UTD-LR Apex')
    xlabel('W_1 (m)')
    ylabel('Mean absolute error (dBA)')
    
    saveas(gcf, [saveDir filesep 'HODComparison_' num2str(k)], 'epsc')
    saveas(gcf, [saveDir filesep 'HODComparison_' num2str(k)], 'svg')

    color = colorStore([1,2,6,3],:);
    figure
    plot(x, meanNN1)
    hold on
    grid on
    plot(x, meanUTD1)
    plot(x, meanNN2, '--')
    plot(x, meanUTD2, '--')
    colororder(color)
    xlim([0.1 3])
    ylim([0 3])
    legend('NN-IIR Extension 1', 'UTD-LR-I Extension 1', 'NN-IIR Extension 2', 'UTD-LR-I Extension 2')
    title(num2str(store(k)))
    xlabel('W_1 (m)')
    ylabel('Mean absolute error (dBA)')
    
    saveas(gcf, [saveDir filesep 'HODComparison_1v2_' num2str(k)], 'epsc')
    saveas(gcf, [saveDir filesep 'HODComparison_1v2_' num2str(k)], 'svg')

    color = colorStore([4,1,2,5,3],:);

    figure
    plot(x, meanBE)
    hold on
    grid on
    plot(x, meanNE)
    plot(x, meanUE)
    plot(x, meanBEI, '--')
    plot(x, meanUEI, '--')
    colororder(color)
    xlim([0.1 3])
    ylim([0 3])
    legend('BTM Extension', 'NN-IIR Extension', 'UTD-LR Extension', 'BTM-I Extension', 'UTD-LR-I Extension')
    title(num2str(store(k)))
    xlabel('W_1 (m)')
    ylabel('Mean absolute error (dBA)')
    
    saveas(gcf, [saveDir filesep 'HODComparison_I_' num2str(k)], 'epsc')
    saveas(gcf, [saveDir filesep 'HODComparison_I_' num2str(k)], 'svg')
end
avErrorN = sum(meanNE) / numBins;
avErrorU = sum(meanUE) / numBins;

[meanBA, meanBE, meanNE, meanNA, meanUE, meanUA, meanUEI] = deal(zeros(1, 2 * numBins));
x = zeros(1, 2 * numBins + 1);
for i = 1:numBins

    idx = width * (i - 1) < W & W <= width * (i);

    i1 = 2 * i - 1;
    i2 = 2 * i;
    meanBA(i1) = mean(meanBtmApex(idx));
    meanBE(i1) = mean(meanBtmExt(idx));
    meanNA(i1) = mean(meanNNApex(idx));
    meanNE(i1) = mean(meanNNExt(idx));
    meanUA(i1) = mean(meanUtdApex(idx));
    meanUE(i1) = mean(meanUtdExt(idx));
    meanUEI(i1) = mean(meanUtdIExt(idx));

    meanBA(i2) = mean(meanBtmApex(idx));
    meanBE(i2) = mean(meanBtmExt(idx));
    meanNA(i2) = mean(meanNNApex(idx));
    meanNE(i2) = mean(meanNNExt(idx));
    meanUA(i2) = mean(meanUtdApex(idx));
    meanUE(i2) = mean(meanUtdExt(idx));
    meanUEI(i2) = mean(meanUtdIExt(idx));
    x(i2:i2 + 1) = width * i;
end
x = x(1:end - 1);
idx = isnan(meanBA);
[meanBA(idx), meanBE(idx), meanNA(idx), meanNE(idx), meanUA(idx), meanUE(idx), meanUEI] = deal(0);

meanBA(idx) = 0;
meanBE(idx) = 0;
meanNA(idx) = 0;
meanNE(idx) = 0;
meanUA(idx) = 0;
meanUE(idx) = 0;
meanUEI(idx) = 0;

figure
plot(x, meanBE)
hold on
grid on
plot(x, meanNE)
plot(x, meanUE)
plot(x, meanBA, '--')
plot(x, meanNA, '--')
plot(x, meanUA, '--')
plot(x, meanUEI)
colororder(color)
xlim([0 3])
legend('BTM Extension', 'NN-IIR Extension', 'UTD-LR Extension', 'BTM Apex', 'NN-IIR Apex', 'UTD-LR Apex', 'UTDI-LR Extension')
xlabel('W_1 (m)')
ylabel('Mean absolute error (dBA)')

saveas(gcf, [saveDir filesep 'HODComparisonSteps'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparisonSteps'], 'svg')


%% Figures
close all

meanBtmApexTot = mean(meanBtmApex);
meanBtmExtTot = mean(meanBtmExt);
meanUtdApexTot = mean(meanUtdApex);
meanUtdExtTot = mean(meanUtdExt);

countBtmExtApex = 100 * sum(meanBtmExt < meanBtmApex) / numPaths;
countBtmUtdApex = 100 * sum(meanBtmExt < meanUtdApex) / numPaths;
countBtmUtdExt = 100 * sum(meanBtmExt < meanUtdExt) / numPaths;

[N, edges, bin] = histcounts(sqrt(meanUtdApex), 50);

edges = linspace(0, 7, 50);
figure
histogram(sqrt(meanBtmApex), length(N), 'BinEdges',edges)
title('BTM Apex')
ylim([0 30])

figure
histogram(sqrt(meanBtmExt), length(N), 'BinEdges',edges)
title('BTM Ext')
ylim([0 30])

figure
histogram(sqrt(meanUtdApex), length(N), 'BinEdges',edges)
title('UTD Apex')
ylim([0 30])

figure
histogram(sqrt(meanUtdExt), length(N), 'BinEdges',edges)
title('UTD Ext')
ylim([0 30])

%% Scores

scoreBtm = sum(meanBtmExt < meanBtmApex);
scoreUtd = sum(meanBtmExt < meanUtdApex);