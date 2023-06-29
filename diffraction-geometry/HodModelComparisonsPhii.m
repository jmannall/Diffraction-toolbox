close all
clear all

set(0, 'DefaultLineLineWidth', 1.5);
colorStore = colororder('default');

%% Data
fs = 48e3;
nfft = 8192;
c = 344;
numEdges = 2;

h = 5;
%h = 2.4;
[zS, zR] = deal(h / 2);
%[zS, zR] = deal(0.5);

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges, 'c', c, 'saveFiles', 3, 'noDirect', true, 'interpolated', false);

%% Generate geometry and controls

numExamples = 500;
numEdges = 2;
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

epsilon = 1e-2;

%% UTD-LR parameters

input = [1; zeros(5, 1)];
windowLength = length(input);
validPath = true(1, numEdges);


%% Models

disp('Start')

% wI = [270 270
%     270 270
%     270 270];
% thetaS = [5; 5; 5];
% thetaR = [250; 250; 250];
% rS = [1; 1; 1];
% rR = [3; 3; 3];
% W = [0.8; 0.8; 0.8];

r = [0.1 3];
rS = RandomUniformDistribution(r, numExamples);
rR = RandomUniformDistribution(r, numExamples);
W = RandomUniformDistribution(r, numExamples);
mA = deg2rad(epsilon) * ones(1, numEdges);


phii = [90; 45; 20];
p = [atand(W / h) 90 * ones(size(W))];
phii = RandomUniformDistribution(p, numExamples);
len = rS + rR + W;
test = tand(phii);
z = h - W ./ tand(phii);
zRange = (h - z) / 2 + [zeros(size(z)) z];
zA = 2.5;
zA = RandomUniformDistribution(zRange, numExamples);
z = tand(90 - phii);
zS = zA - (W / 2 + rS) .* z;
zR = zA + (W / 2 + rR) .* z;

dZ = zS - zR;
L = sqrt((rS + sum(W, 2) + rR) .^ 2 + dZ .^ 2);
%%
[zA1, phii1] = CalculateApex(rS', rR' + W', zS', zR', h, true);
[zA2, phii2] = CalculateApex(rS' + W', rR', zS', zR', h, true);

zAall = [zA1; zA2];
if max(zAall(:,3)) > h || min(zAall(:,3)) < 0
    disp('Apex outside of edge')
    return
end
createPlot = true;
n = 1;
%%

index = '2180cfff4ddad1de505214480a4535ce';
[loadPath, savePath] = deal(['geometry/NthOrderPaths_01mTo3m_', num2str(index), '.mat']);
load(loadPath, 'data')

%%
n = 1;
load('HODdata.mat')
zA1 = zA1(:,3);
zA2 = zA2(:,3);


numExamples = 500;
for i = n:numExamples

%     wI(i,:) = 185 + (355 - 185) * rand(1, numEdges);
%     thetaS(i) = 0.1 + (wI(i,1) - 180) * rand(1);
%     thetaR(i) = 180 + (wI(i,end) - 180.1) * rand(1);
%     bA(i,:) = deg2rad([wI(i,1) - thetaS(i), wI(i,2:numEdges - 1), thetaR(i)]);

    wI(i,:) = data(i).wedgeIndex;
    thetaS(i) = data(i).thetaS;
    thetaR(i) = data(i).thetaR;
    bA(i,:) = deg2rad([wI(i,1) - thetaS(i), wI(i,2:numEdges - 1), thetaR(i)]);

    [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid, vReceiver] = CreateNthOrderPathData(wI(i,:), thetaS(i), thetaR(i), rS(i,:), rR(i,:), W(i,:), h);

    source(:,3) = zS(i);
    receiver(:,3) = zR(i);
    controlparameters.fs = 2 * fs;
    controlparameters.nfft = 2 * nfft;
    [~, tfmagStore, ~, fvecBtm, ~] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);
    tfmag.Btm = tfmagStore.diff2;
    tfmagN.Btm(:,i) = CreateNBandMagnitude(tfmag.Btm, fidx);

    if zS(i) < zR(i)
        z = zS(i);
    else
        z = zR(i);
    end
    dZ = zR(i) - zS(i);

    apex(1,3) = zS(i) + (rS(i)) * dZ / (rS(i) + W(i) + rR(i));
    apex(2,3) = zS(i) + (rS(i) + W(i)) * dZ / (rS(i) + W(i) + rR(i));

    if max(apex(:,3)) > h || min(apex(:,3)) < 0
        disp('Apex outside of edge')
    end

    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    rSE = [rS(i), W(i)];
    rRE = [rR(i), fliplr(W(i))];

    cumRs = cumsum(rSE)';
    cumRr = fliplr(cumsum(rRE))';

    vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
    vector = vSource - vCorners;
    vector(:,3) = 0;
    vector = cumRs .* vector ./ vecnorm(vector, 2, 2);
    vSource = vCorners + vector;
    vSource(:,3) = zS(i);

    vector = vReceiver - vCorners;
    vector(:,3) = 0;
    vector = cumRr .* vector ./ vecnorm(vector, 2, 2);
    vReceiver = vCorners + vector;
    vReceiver(:,3) = zR(i);

    dZ = cumsum(abs([zS(i) - zA1(i); zA1(i) - zA2(i)]));
    dS = sqrt(cumRs .^ 2 + dZ .^ 2);
    dZ = flipud(cumsum(abs([zA2(i) - zR(i); zA1(i) - zA2(i)])));
    dR = sqrt(cumRr .^ 2 + dZ .^ 2);
    [tfmag.BtmE, ~, tf.BtmE] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, dS, dR, L(i));

    tfmagN.BtmE(:,i) = CreateNBandMagnitude(tfmag.BtmE(:,end), fidx);

    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    dSA = [sqrt(rS(i) ^ 2 + (zS(i) - zA1(i)) ^ 2), sqrt(W(i) ^ 2 + (zA1(i) - zA2(i)) ^ 2)];
    dRA = [sqrt(W(i) ^ 2 + (zA1(i) - zA2(i)) ^ 2), sqrt(rR(i) ^ 2 + (zR(i) - zA2(i)) ^ 2)];

    [tfmag.BtmA, ~, tf.BtmA] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, dSA, dRA, L(i));

    tfmagN.BtmA(:,i) = CreateNBandMagnitude(tfmag.BtmA(:,end), fidx);

    %% Interpolated

    [tfmag.BtmIE, tfmag.BtmIA, tf.BtmIE, tf.BtmIA] = deal(zeros(nfft, numEdges + 1));

    disp('BTM Ext')

    [tfmag.BtmIE(:,1), ~, tf.BtmIE(:,1)] = SingleWedgeInterpolated(h, wI(i,1), epsilon, wI(i,1) - thetaS(i), rS(i), W(i) + rR(i), zS(i), zR(i), controlparameters, false);
    [tfmag.BtmIE(:,2), ~, tf.BtmIE(:,2)] = SingleWedgeInterpolated(h, wI(i,2), epsilon, thetaR(i), rS(i) + W(i), rR(i), zS(i), zR(i), controlparameters, false);
    tf.BtmIE(:,end) = 0.5^(numEdges - 1) * (1 / L(i)) * prod(tf.BtmIE(:,1:numEdges), 2);
    tfmag.BtmIE(:,end) = mag2db(abs(tf.BtmIE(:,end)));

    tfmagN.BtmIE(:,i) = CreateNBandMagnitude(tfmag.BtmIE(:,end), fidx);

    disp('BTM Apex')

    [tfmag.BtmIA(:,1), ~, tf.BtmIA(:,1)] = SingleWedgeInterpolated(h, wI(i,1), epsilon, wI(i,1) - thetaS(i), rS(i), W(i), zS(i), zA2(i), controlparameters, false);
    [tfmag.BtmIA(:,2), ~, tf.BtmIA(:,2)] = SingleWedgeInterpolated(h, wI(i,2), epsilon, thetaR(i), W(i), rR(i), zA1(i), zR(i), controlparameters, false);
    tf.BtmIA(:,end) = 0.5^(numEdges - 1) * (1 / L(i)) * prod(tf.BtmIA(:,1:numEdges), 2);
    tfmag.BtmIA(:,end) = mag2db(abs(tf.BtmIA(:,end)));

    tfmagN.BtmIA(:,i) = CreateNBandMagnitude(tfmag.BtmIA(:,end), fidx);

    %% NN

    controlparameters.fs = fs;
    controlparameters.nfft = nfft;
    [tfmag.NNE, tfmag.NNA, tf.NNE, tf.NNA] = deal(zeros(nfft / 2, numEdges + 1));

    disp('NN Ext')

    r = [rS(i), rS(i) + W(i); W(i) + rR(i), rR(i)];
    r1 = min(r);
    r2 = max(r);

    [~,idx] = min(r);
    [z1, z2] = deal(zeros(1,2));
    for j = 1:numEdges
        if idx(j) == 1
            z1(j) = zS(i);
            z2(j) = zR(i);
        else
            z1(j) = zR(i);
            z2(j) = zS(i);
        end
        if z1(j) > h / 2
            z1(j) = h - z1(j);
            z2(j) = h - z2(j);
        end
    end


    X = [deg2rad(wI(i,:)); bA(i,:); mA; [h, h]; r1; r2; z1; z2];
    %X = [deg2rad(wI); bA; mA; wL / 2; r1; r2; [0.01, 0.01]; [0.01, 0.01]];
    X = dlarray(single(X), "CB");

    [tfmag.NNE(:,1:numEdges), ~, tf.NNE(:,1:numEdges)] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);

    tf.NNE(:,end) = 0.5^(numEdges - 1) * (1 / L(i)) * prod(tf.NNE(:,1:numEdges), 2);
    tfmag.NNE(:,end) = mag2db(abs(tf.NNE(:,end)));

    tfmagN.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,end), fidx);

    disp('NN Apex')

    r = [rS(i), W(i); W(i), rR(i)];
    r1 = min(r);
    r2 = max(r);

    [~,idx] = min(r);
    [z1, z2] = deal(zeros(1,2));

    if idx(1) == 1
        z1(1) = zS(i);
        z2(1) = zA2(i);
    else
        z1(1) = zA2(i);
        z2(1) = zS(i);
    end
    if idx(2) == 1
        z1(2) = zA1(i);
        z2(2) = zR(i);
    else
        z1(2) = zR(i);
        z2(2) = zA1(i);
    end
    for j = 1:numEdges
        if z1(j) > h / 2
            z1(j) = h - z1(j);
            z2(j) = h - z2(j);
        end
    end

    X = [deg2rad(wI(i,:)); bA(i,:); mA; [h, h]; r1; r2; z1; z2];
    %X = [deg2rad(wI); bA; mA; wL / 2; r1; r2; [0.01, 0.01]; [0.01, 0.01]];
    X = dlarray(single(X), "CB");

    [tfmag.NNA(:,1:numEdges), ~, tf.NNA(:,1:numEdges)] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);

    tf.NNA(:,end) = 0.5^(numEdges - 1) * (1 / L(i)) * prod(tf.NNA(:,1:numEdges), 2);
    tfmag.NNA(:,end) = mag2db(abs(tf.NNA(:,end)));

    tfmagN.NNA(:,i) = CreateNBandMagnitude(tfmag.NNA(:,end), fidx);

    %% UTD

    disp('UTD Ext')

    tfcomplex = zeros(4, numEdges + 1);

    for j = 1:numEdges
        [~, ~, tfcomplexStore] = SingleUTDWedgeInterpolated(epsilon, rad2deg(mA(j) + bA(i,j)), cumRs(j), cumRr(j), wI(i,j), phii1(i), controlparameters);
        tfcomplex(:,j) = L(i) * tfcomplexStore;
    end
    
    tfcomplex(:,end) = 0.5^(numEdges - 1) .* (1 / L(i)) .* prod(tfcomplex(:,1:numEdges), 2);
    tfmagStore = mag2db(abs(tfcomplex));

    [~, ~, tfmag.UtdILRE] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tfmagStore', c, fs, validPath);
    [tfmagN.UtdILRE(:,i), ~, ~] = CreateFrequencyNBands(tfmag.UtdILRE(:,end), fvec, nBand);
    
    %%
    disp('UTD Apex')

    tfcomplex = zeros(4, numEdges + 1);

    rSA = [rS(i); W(i)];
    rRA = [W(i); rR(i)];
    for j = 1:numEdges
        [~, ~, tfcomplexStore] = SingleUTDWedgeInterpolated(epsilon, rad2deg(mA(j) + bA(i,j)), rSA(j), rRA(j), wI(i,j), phii1(i), controlparameters);
        tfcomplex(:,j) = (dSA(j) + dRA(j)) * tfcomplexStore;
    end
    
    tfcomplex(:,end) = 0.5^(numEdges - 1) .* (1 / L(i)) .* prod(tfcomplex(:,1:numEdges), 2);
    tfmagStore = mag2db(abs(tfcomplex));

    [~, ~, tfmag.UtdILRA] = DelayLineLR(input, zeros(size(pathLength)), windowLength, validPath, tfmagStore', c, fs, validPath);
    [tfmagN.UtdILRA(:,i), ~, ~] = CreateFrequencyNBands(tfmag.UtdILRA(:,end), fvec, nBand);
    
    %% Plot
    close all

%     figure
%     semilogx(fvecBtm, tfmag.Btm)
%     hold on
%     grid on
%     semilogx(fvecBtm, tfmag.BtmE(:,end))
%     semilogx(fvecBtm, tfmag.BtmIE(:,end))
%     semilogx(fvec, tfmag.NNE(:,end))
%     semilogx(fvec, tfmag.NNA(:,end))
%     semilogx(fvec, tfmag.UtdILRE(:,end))
%     semilogx(fvec, tfmag.UtdILRA(:,end))
%     legend('BTM', 'BTME', 'BTMIE', 'NNE', 'NNA', 'UtdILR', 'UtdILRA')
%     xlim([20 20e3])

    tfmagN1.BtmIE(:,i) = CreateNBandMagnitude(tfmag.BtmIE(:,1), fidx);
    tfmagN1.UtdILRE(:,i) = CreateFrequencyNBands(tfmag.UtdILRE(:,1), fvec, nBand);
    tfmagN1.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,1), fidx);

    tfmagN2.BtmIE(:,i) = CreateNBandMagnitude(tfmag.BtmIE(:,2), fidx);
    tfmagN2.UtdILRE(:,i) = CreateFrequencyNBands(tfmag.UtdILRE(:,2), fvec, nBand);
    tfmagN2.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,2), fidx);

    n = n + 1;
    disp(num2str(i))
end

%%

loss = CalculateLoss(tfmagN, tfmagN.Btm);
loss1 = CalculateLoss(tfmagN1, tfmagN1.BtmIE);
loss2 = CalculateLoss(tfmagN2, tfmagN2.BtmIE);

%% Plot

close all

[test, testIdx] = sort(loss.i.BtmIA);

if isfield(loss, 'w')
    loss = rmfield(loss, 'w');
end

test = sqrt((zA1 - zA2) .^ 2 + W .^ 2);

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
    num(i) = sum(idx);
end
for j = 1:numFields
    field = fields{j};
    idx = isnan(loss.w.(field));
    loss.w.(field)(idx) = 0;
    loss.w.(field)(numBins + 2) = 0;
end

saveDir = 'figures';
color = colorStore([4,4,1,2,4,4,1,2],:);

figure
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
ylabel('Mean dB difference (dB)')

saveas(gcf, [saveDir filesep 'HODComparisonW'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparisonW'], 'svg')

%% Phii

if isfield(loss, 'p')
    loss = rmfield(loss, 'p');
end

test = sqrt((zA1 - zA2) .^ 2 + W .^ 2);

fields = fieldnames(loss.i);
numFields = length(fields);

width = 3;
numBins = 90 / width;
x = (0:width:90 + width) - width / 2;
for i = 2:numBins+1
    idx = width * (i - 2) < phii & phii <= width * (i - 1);
    for j = 1:numFields
        field = fields{j};
        loss.p.(field)(i) = mean(loss.i.(field)(idx));
    end
    num(i) = sum(idx);
end
for j = 1:numFields
    field = fields{j};
    idx = isnan(loss.p.(field));
    loss.p.(field)(idx) = 0;
    loss.p.(field)(numBins + 2) = 0;
end

saveDir = 'figures';
color = colorStore([4,4,1,2,4,4,1,2],:);

figure
plot(x, loss.p.BtmE, 'Color', [color(1,:), 0.6])
hold on
grid on
plot(x, loss.p.BtmIE)
plot(x, loss.p.NNE)
plot(x, loss.p.UtdILRE)
plot(x, loss.p.BtmA, '--', 'Color', [color(1,:), 0.6])
plot(x, loss.p.BtmIA, '--')
plot(x, loss.p.NNA, '--')
plot(x, loss.p.UtdILRA, '--')
colororder(color)
xlim([1 90])
ylim([0 15])
legend('BTM Extension', 'BTM-I Extension', 'NN-IIR (best) Extension', 'UTD-LR Extension', 'BTM Apex', 'BTM-I Apex', 'NN-IIR (best) Apex', 'UTD-LR Apex')
xlabel('phii (degrees)')
ylabel('Mean dB difference (dB)')

saveas(gcf, [saveDir filesep 'HODComparisonPhii'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparisonPhii'], 'svg')

%% D

if isfield(loss, 'd')
    loss = rmfield(loss, 'd');
end

D = sqrt((zA1 - zA2) .^ 2 + W .^ 2);

fields = fieldnames(loss.i);
numFields = length(fields);

width = 0.2;
numBins = 6 / width;
x = (0:width:6 + width) - width / 2;
for i = 2:numBins+1
    idx = width * (i - 2) < D & D <= width * (i - 1);
    for j = 1:numFields
        field = fields{j};
        loss.d.(field)(i) = mean(loss.i.(field)(idx));
    end
    num(i) = sum(idx);
end
for j = 1:numFields
    field = fields{j};
    idx = isnan(loss.w.(field));
    loss.d.(field)(idx) = 0;
    loss.d.(field)(numBins + 2) = 0;
end

saveDir = 'figures';
color = colorStore([4,4,1,2,4,4,1,2],:);

figure
plot(x, loss.d.BtmE, 'Color', [color(1,:), 0.6])
hold on
grid on
plot(x, loss.d.BtmIE)
plot(x, loss.d.NNE)
plot(x, loss.d.UtdILRE)
plot(x, loss.d.BtmA, '--', 'Color', [color(1,:), 0.6])
plot(x, loss.d.BtmIA, '--')
plot(x, loss.d.NNA, '--')
plot(x, loss.d.UtdILRA, '--')
colororder(color)
xlim([0.1 5])
ylim([0 6])
legend('BTM Extension', 'BTM-I Extension', 'NN-IIR (best) Extension', 'UTD-LR Extension', 'BTM Apex', 'BTM-I Apex', 'NN-IIR (best) Apex', 'UTD-LR Apex', 'Location', 'northwest')
xlabel('W_1 (m)')
ylabel('Mean dB difference (dB)')

saveas(gcf, [saveDir filesep 'HODComparisonD'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparisonD'], 'svg')