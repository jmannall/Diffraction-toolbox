
clear all
close all
set(0, 'DefaultLineLineWidth', 1.5);

fvec = 1000;
c = 344;
k = 2 * pi * fvec/ c; % Precompute
pathLength = 2;
DirRef1 = mag2db(abs(exp(-1i * k * pathLength) / pathLength));
z = 0.1 + 0.7i;
test1 = mag2db(abs(z));
test2 = 20 * log10(sqrt(real(z) ^ 2 + imag(z) ^ 2));
DirRef2 = mag2db(abs(1 / pathLength));

%% Terms

%{

Path finding cost (Cp): Cost of finding all the relevant paths every update
Calculation cost (Cc): Cost of calculating UTD / infering from a NN / 
calculating coefficients for a filter structure every update

Rendering cost (Cr): Cost of running a filter structure in FLOPS

Diffraction model cost (Cd) = Fr * Cc + Cr
GA model cost (Cga) = Fr * (Cp + Cc) + Fs * Cr

Update rate (Fr)
Sampling frequency (Fs)

%}

%% Constants

Fs = 48e3;
Fr = 0:100;
savePath = 'figures/EURASIP';

colorStore = colororder;

%% Filters cost per sample

opPerModule = 17;
numModules = 10;
numSum = 5;
numMult = 4;
LR = opPerModule * numModules + numSum + numMult;
numSum = 4;
numMult = 5;
IIR2 = numSum + numMult;
numSum = 8;
numMult = 9;
IIR4 = numSum + numMult;

%% Filtering running costs

diffOrders = 1:24;
numDiff = length(diffOrders);
CrLR = Fs * LR * ones(1, numDiff);
CrIIR2 = Fs * (IIR2 + (IIR2 - 1) * diffOrders);
CrIIR4 = Fs * (IIR4 + (IIR4 - 1) * diffOrders);

figure
plot(diffOrders, CrIIR2)
hold on
plot(diffOrders, CrLR)
plot(diffOrders, CrIIR4)
grid on
legend('2nd order IIR', 'L-R', '4th order IIR', 'Location', 'southeast')
xlabel('Diffraction order')
ylabel('Rendering cost (MFLOPS)')
yticklabels(0:10)
title({'Rendering cost of each filter structure', 'against diffraction order'})

saveas(gcf, [savePath, filesep, 'RenderingCosts.svg'], 'svg')

diffOrders = [1 3 5 7 10];
numDiff = length(diffOrders);

CrLR = Fs * LR * ones(1, numDiff);
CrIIR2 = Fs * (IIR2 + (IIR2 - 1) * diffOrders);
CrIIR4 = Fs * (IIR4 + (IIR4 - 1) * diffOrders);

%% UTD calculation costs

numUtdFreq = 4;
utdPerFreq = 58;
utdOnce = 66;
tfmagCalc = 8;
dirRef = 4 + tfmagCalc;
interpolate = 6 + numUtdFreq * 5;

CcUtd = 2 * (utdOnce + numUtdFreq * utdPerFreq) + dirRef + interpolate;
CcLR = numUtdFreq;

CcUtdLR = CcUtd + CcLR;

%% NN calculation costs

gx = 5; % Activation layer and normalisation cost per node
nN = 25; % Number of nodes per layer
nL = 3; % Number of layers
nO = 5; % Number of outputs
nI = 8; % Number of inputs

gx = 5; % Activation layer and normalisation cost per node
nN = 36; % Number of nodes per layer
nL = 7; % Number of layers
nO = 5; % Number of outputs
nI = 8; % Number of inputs

CcPZTransform = 6;
numPZ = 4;
CcKTransform = 1;
CcOutputTransform = numPZ * CcPZTransform + CcKTransform;

nNTotal = nN * nL + nO; % Total number of nodes in DNN (i.e the number of bias additions)
nNActivations = nN * nL; % Total number of nodes in DNN that require activation/normalisation
nCTotal = nI * nN + nN ^ 2 * (nL - 1) + nN * nO; % Total number of connections in DNN (i.e the number of weights)

CcNN = nNTotal + nNActivations * gx + nCTotal + CcOutputTransform;
CcIIR = 6;

CcNNIIR = CcNN + CcIIR;

gx = 5; % Activation layer and normalisation cost per node
nN = 20; % Number of nodes per layer
nL = 4; % Number of layers
nO = 5; % Number of outputs
nI = 8; % Number of inputs

CcOutputTransform = numPZ * CcPZTransform + CcKTransform;

nNTotal = nN * nL + nO; % Total number of nodes in DNN (i.e the number of bias additions)
nNActivations = nN * nL; % Total number of nodes in DNN that require activation/normalisation
nCTotal = nI * nN + nN ^ 2 * (nL - 1) + nN * nO; % Total number of connections in DNN (i.e the number of weights)

CcNN2 = nNTotal + nNActivations * gx + nCTotal + CcOutputTransform;

CcNNIIR2 = CcNN2 + CcIIR;

%% UDFA calculation costs

n = 4;
numUDFAt = n + 1;
numUDFAi = n;
UDFAOnce = 88;
UDFAPerFreq = 10;

CcUDFA = UDFAOnce + UDFAPerFreq * (numUDFAt + numUDFAi);

CcFsh = 15; % Fsh adjustment
CcFc = CcFsh * numUDFAi;
CcIIR = 12; % Calculate 1 filters a b parameters
CcIIRs = (CcIIR + 1) * n + 1; % prod(b1_n) + (mult gt(1) / L) to first gain parameter
CcUDFAIIR = CcUDFA + CcFc + CcIIRs;

%% BTM calculation costs

single_thread_flops = 3.2*10^9;
numPaths = 140;
%tb = 3.31*10^-3;
tb = 137.51*10^-3;

CcBtmOld = tb * single_thread_flops / numPaths;

wedgeLength = 5;
r1 = 1;
r2 = 1;
z1 = wedgeLength / 2;
z2 = wedgeLength / 2;
t0 = (r1 + r2) / c;
tMax = (sqrt(r1 ^ 2 + z1 ^ 2) + sqrt(r2 ^ 2 + z2 ^ 2)) / c;
numSegments = 2 * ceil((tMax - t0) * Fs);
numTheta = 4;
numIntegralPoints = 2 * numSegments + 1;

btmOnce = 3;
btmPerTheta = 3;
btmEveryIntegralPoint = 9;
btmEveryIntegralPointTheta = 3;
btmEverySegment = 4;
btmEverySegementTheta = 8;
CcSimpInt = 8; % 3-point simpsons rule

CcBtm = btmOnce + btmPerTheta * numTheta + btmEveryIntegralPoint * numIntegralPoints + btmEveryIntegralPointTheta * numIntegralPoints * numTheta + btmEverySegment * numSegments + CcSimpInt * numSegments * numTheta;

% Calculating subdivisions
% n = fs * (sqrt(r1 ^ 2 + z1 ^ 2) + sqrt(r2 ^ 2 + z2 ^ 2) - (r1 + r2)) / c;
% n = numEdgeDivisions - 1;
% n * c / fs + (r1 + r2) = sqrt(r1 ^ 2 + (z - z1) ^ 2) + sqrt(r2 ^ 2 + (z - z2) ^ 2)


%% Running costs of each method

CdUtdLR = Fr' * diffOrders .* CcUtdLR + CrLR;
CdNNIIR = Fr' * diffOrders .* CcNNIIR + CrIIR2;
CdNNIIR2 = Fr' * diffOrders .* CcNNIIR2 + CrIIR2;

%CrIIR4 = ones(size(CrIIR4)) * CrIIR4(1);
CdUDFAIIR = Fr' * diffOrders .* CcUDFAIIR + CrIIR4;

color = colorStore(1:5,:);
figure
colororder(color);
plot(Fr, CdNNIIR)
hold on
plot(Fr, CdUtdLR, '--')
plot(Fr, CdUDFAIIR, '-.')
legend('NN-IIR', 'UTD-LR', 'UDFA-IIR')
xlabel('Update rate (Hz)')
ylabel('Cost (MFLOPS)')
l = legend(string(diffOrders), 'Location','northwest');
title(l, 'Diffraction Order')
%yticklabels(0:2:14)
title('Diffraction model cost')

%saveas(gcf, [savePath, filesep, 'ModelCosts.svg'], 'svg')

%% Overlap add methods

T60 = 0.03;
N = Fs./Fr'+T60*Fs-1;

numUtdFreq = N;
numUtdFreq = 512 * ones(size(N));
utdPerFreq = 84;
utdOnce = 67;
tfmagCalc = 8;
dirRef = 4 + tfmagCalc;
interpolate = 6 + numUtdFreq * 5;

CcUtd = 2 * (utdOnce + numUtdFreq * utdPerFreq) + dirRef + interpolate;
CcUtd = utdOnce + numUtdFreq * utdPerFreq;

CrOverlapAddUtd = Fr'.*(12*N.*log2(N)+6*N+T60*Fs-1);
CrOverlapAddBtm = Fr'.*(18*N.*log2(N)+6*N+T60*Fs-1);

CdNNIIR = Fr' * diffOrders .* CcNNIIR + CrIIR2;

CdUtdOverlapAdd = Fr' * diffOrders .* CcUtd + CrOverlapAddUtd;
CdBtmOverlapAdd = Fr' .* CcBtm + CrOverlapAddBtm;

%% Diffraction model costs

color = colorStore(1:5,:);
figure
colororder(color);
semilogy(Fr, CdNNIIR)
hold on
grid on
semilogy(Fr, CdUtdLR, '--')
semilogy(Fr, CdUDFAIIR, '-.')
semilogy(Fr, CdUtdOverlapAdd, '-.')
semilogy(Fr, CdBtmOverlapAdd, ':')
legend('NN-IIR', 'UTD-LR', 'UDFA-IIR', 'UTD-Overlap add', 'BTM-Overlap add')
xlabel('Update rate (Hz)')
ylabel('Cost (FLOPS)')
l = legend(string(diffOrders), 'Location','northwest');
title(l, 'Diffraction Order')
%yticklabels(0:2:14)
ylim([5e5 2e9])
xlim([10 100])
title('Diffraction model cost')

%% Path finding complexity

close all

single_thread_flops = 3.5*10^9;

% 6 examples taken from table 1 Schissler 2021
% frameTime = [2.1 36.5 0.28 4.2 120.7 72] * 10^-3; % [s]
% numSources = [6 7 2 6 16 20];

% 6 scenes taken from figure 10 Schissler 2021
appartment = [0.178, 0.231, 0.213, 0.248, 0.303] * 10^-3;
sibenik = [0.118, 0.2, 0.353, 0.905, 0.982] * 10^-3;
sunTemple = [1.303, 2.701, 3.899, 4.875, 6.392] * 10^-3;
bistro = [1.466, 2.459, 3.634, 4.651, 6.098] * 10^-3;
sponza = [0.127, 0.256, 0.321, 0.344, 0.361] * 10^-3;
warehouse = [0.814, 1.831, 2.974, 3.762, 4.762] * 10^-3;

frameTime = [appartment; sibenik; sunTemple; bistro; sponza; warehouse];
% 3 examples taken from table 1 Schissler 2021
% frameTime = [2.1 0.28 4.2] * 10^-3; % [s]
% numSources = [6 2 6];

availableProcessing = single_thread_flops - CdUtdLR;

numScenes = size(frameTime, 1);
Cp = zeros(length(Fr), numScenes, numDiff);
for i = 1:numScenes
    oneUpdate = single_thread_flops .* frameTime(i,:); % FLOPS required to calculate one update for each scene
    Cp(:,i,:) = Fr' * oneUpdate;
    % Accounting for processing taken up by UtdLR calculations
    oneUpdate = availableProcessing .* frameTime(i,:); % FLOPS required to calculate one update for each scene
    CpAP(:,i,:) = Fr' .* oneUpdate;
end

x = 300;
gap = 10;

color = colorStore(1:6,:);
figure('Position', [gap gap 6 * x 3 * x])
t = tiledlayout(2, 3, 'TileSpacing', 'compact');
title(t, 'Path finding and diffraction model costs')
for i = 1:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, Cp(:,:,i))
    hold on
    semilogy(Fr, CdUtdLR(:,i), '--')
    semilogy(Fr, CdNNIIR(:,i), '--')
    semilogy(Fr, CdUDFAIIR(:,i), '--')
    semilogy(Fr, CdUtdOverlapAdd(:,i), '--')
    if i == 1
        semilogy(Fr, CdBtmOverlapAdd(:,i), '--')
    end
    l = legend('apartment', 'sibenik', 'sun temple', 'bistro', 'sponza', 'warehouse', 'UTD-LR', 'NN-IIR', 'UTD-Overlap add', 'BTM-Overlap add', 'Location', 'southeast');
    title(l, 'Scene and diffraction model')
    xlim([10 100])
    ylim([5e5 3e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end

%saveas(gcf, [savePath, filesep, 'ComponentCosts.svg'], 'svg')

figure('Position', [gap gap 6 * x 3 * x])
t = tiledlayout(2, 3, 'TileSpacing', 'compact');
title(t, 'Path finding costs')
for i = 1:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, Cp(:,:,i))
    l = legend('apartment', 'sibenik', 'sun temple', 'bistro', 'sponza', 'warehouse', 'Location', 'southeast');
    title(l, 'Scene')
    xlim([10 100])
    ylim([1e6 3e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end

%saveas(gcf, [savePath, filesep, 'PathFindingCosts.svg'], 'svg')

figure('Position', [gap gap 6 * x 3 * x])
t = tiledlayout(2, 3, 'TileSpacing', 'compact');
title(t, 'GA model costs')
for i = 1:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, Cp(:,:,i) + CdNNIIR(:,i))
    hold on
    semilogy(Fr, Cp(:,:,i) + CdUtdLR(:,i), '--')
    semilogy(Fr, Cp(:,:,i) + CdUDFAIIR(:,i), '--')
    semilogy(Fr, Cp(:,:,i) + CdUtdOverlapAdd(:,i), '-.')
    if i == 1
        semilogy(Fr, Cp(:,:,i) + CdBtmOverlapAdd(:,i), ':')
    end
    l = legend('apartment', 'sibenik', 'sun temple', 'bistro', 'sponza', 'warehouse', 'Location', 'southeast');
    title(l, 'Scene')
    xlim([10 100])
    ylim([1e6 3e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end

%saveas(gcf, [savePath, filesep, 'GAModelCosts.svg'], 'svg')

%% Small plots

close all

textScale = 1.2;

color = colorStore(1:6,:);
figure('Position', [gap gap 6 * x 1.5 * x])
tiledlayout(1, 3, 'TileSpacing', 'compact');
for i = 1:2:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, Cp(:,:,i))
    grid on
    xlim([10 100])
    ylim([1e6 3e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end
l = legend('Apartment', 'Sibenik', 'Sun temple', 'Bistro', 'Sponza', 'Warehouse', 'Location', 'southeast');
title(l, 'Scene')
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, 'PathFindingCostsSmall'], 'epsc')
saveas(gcf, [savePath, filesep, 'PathFindingCostsSmall'], 'svg')

idx = [1 4];
color = colorStore(1:length(idx),:);
figure('Position', [gap gap 6 * x 1.5 * x])
t = tiledlayout(1, 3, 'TileSpacing', 'compact');
for i = 1:2:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, Cp(:,idx,i) + CdNNIIR(:,i))
    hold on
    grid on
    semilogy(Fr, Cp(:,idx,i) + CdUtdLR(:,i), '--')
    semilogy(Fr, Cp(:,idx,i) + CdUDFAIIR(:,i), '-.')
    semilogy(Fr, Cp(:,idx,i) + CdUtdOverlapAdd(:,i), '-.')
    if i == 1
        semilogy(Fr, Cp(:,idx,i) + CdBtmOverlapAdd(:,i), ':')
    end
    xlim([10 100])
    ylim([5e6 5e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end
fontsize(gcf,scale=textScale)
l = legend('Apartment', 'Bistro', 'Location', 'southeast');
title(l, 'Scene')
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, 'GAModelCostsSmall'], 'epsc')
saveas(gcf, [savePath, filesep, 'GAModelCostsSmall'], 'svg')

%%

color = colorStore([1, 1, 2, 5, 3, 4],:);
figure('Position', [gap gap 6 * x 1.5 * x])
t = tiledlayout(1, 3, 'TileSpacing', 'compact');
for i = 1:2:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, CdNNIIR(:,i))
    hold on
    grid on
    semilogy(Fr, CdNNIIR2(:,i), 'Color', [color(1,:), 0.6])
    semilogy(Fr, CdUtdLR(:,i), '-.')
    semilogy(Fr, CdUDFAIIR(:,i), ':')
    semilogy(Fr, CdUtdOverlapAdd(:,i), '--')
    if i == 1
        semilogy(Fr, CdBtmOverlapAdd, ':')
    else
        semilogy(Fr, 1e3 * CdBtmOverlapAdd, ':')
    end
    %semilogy(Fr, CdUtdIIR(:,i), ':')
    xlim([10 100])
    ylim([5e5 5e8])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end
fontsize(gcf,scale=textScale)
l = legend('NN-IIR (best)', 'NN-IIR (small)', 'UTD-LR', 'UDFA', 'UTD-overlap add', 'BTM-overlap add', 'Location', 'bestoutside');
title(l, 'Diffraction Model')
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, 'DiffractionModelCostsSmall'], 'epsc')
saveas(gcf, [savePath, filesep, 'DiffractionModelCostsSmall'], 'svg')

%%
close all

color = colorStore([1, 1, 2, 5, 4],:);
%figure('Position', [gap gap 6 * x 1.5 * x])
figure
colororder(color);
semilogy(Fr, CdNNIIR(:,1))
hold on
grid on
semilogy(Fr, CdNNIIR2(:,1), 'Color', [color(1,:), 0.6])
semilogy(Fr, CdUtdLR(:,1), '--')
semilogy(Fr, CdUDFAIIR(:,1), ':')
%semilogy(Fr, CdUtdOverlapAdd(:,1), '--')
semilogy(Fr, CdBtmOverlapAdd, 'LineWidth', 3)
%semilogy(Fr, CdUtdIIR(:,i), ':')
xlim([10 100])
ylim([5e5 5e8])
xlabel('Update rate (Hz)')
ylabel('FLOPS')
%title(['Diffraction order: ', num2str(diffOrders(i))])
%fontsize(gcf,scale=textScale)
l = legend('NN-IIR (best)', 'NN-IIR (small)', 'UTD-LR', 'UDFA', 'BTM-overlap add', 'Location', 'northwest');
%title(l, 'Diffraction Model')
fontsize(gcf,20,"pixels")

saveas(gcf, [savePath, filesep, '1stOrderDiffractionModelCosts'], 'epsc')
saveas(gcf, [savePath, filesep, '1stOrderDiffractionModelCosts'], 'svg')

return
%%

perNNIIR = reshape(CdNNIIR, 101,1,[]) ./ (Cp + reshape(CdNNIIR, 101,1,[])) * 100;
perUtdLR = reshape(CdUtdLR, 101,1,[]) ./ (Cp + reshape(CdUtdLR, 101,1,[])) * 100;

perUtdOverlapAdd = reshape(CdUtdOverlapAdd, 101,1,[]) ./ (Cp + reshape(CdUtdOverlapAdd, 101,1,[])) * 100;
perBtmOverlapAdd = reshape(CdBtmOverlapAdd, 101,1,[]) ./ (Cp + reshape(CdBtmOverlapAdd, 101,1,[])) * 100;

idx = [1 4];
color = colorStore(1:length(idx),:);
figure('Position', [gap gap 6 * x 1.5 * x])
t = tiledlayout(1, 3, 'TileSpacing', 'compact');
for i = 1:2:numDiff
    nexttile
    colororder(color);
    plot(Fr, perNNIIR(:,idx,i))
    hold on
    grid on
    plot(Fr, perUtdLR(:,idx,i), '--')
    plot(Fr, perUtdOverlapAdd(:,idx,i), '-.')
    if i == 1
        plot(Fr, perBtmOverlapAdd(:,idx,i), ':')
    end
    xlim([10 100])
    %ylim([5e6 5e9])
    %xlabel('Update rate (Hz)')
    %ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end
fontsize(gcf,scale=textScale)
l = legend('Apartment', 'Bistro', 'Location', 'southeast');
title(l, 'Scene')
fontsize(gcf,20,"pixels")

%% Combine plot

close all

idx = [1 4];
color = colorStore(1:length(idx),:);
figure('Position', [gap gap 6 * x 3 * x])
t = tiledlayout(2, 3, 'TileSpacing', 'compact');
for i = 1:2:numDiff
    nexttile
    colororder(color);
    semilogy(Fr, Cp(:,idx,i) + CdNNIIR(:,i))
    hold on
    semilogy(Fr, Cp(:,idx,i) + CdUtdLR(:,i), '--')
    semilogy(Fr, Cp(:,idx,i) + CdUtdOverlapAdd(:,i), '-.')
    if i == 1
        semilogy(Fr, Cp(:,idx,i) + CdBtmOverlapAdd(:,i), ':')
    end
    l = legend('Apartment', 'Bistro', 'Location', 'southeast');
    title(l, 'Scene')
    xlim([10 100])
    ylim([1e6 3e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
end

color = colorStore(1:4,:);
for i = 1:2:numDiff
    nexttile
    semilogy(Fr, CdNNIIR(:,i))
    hold on
    semilogy(Fr, CdUtdLR(:,i), '--')
    semilogy(Fr, CdUtdOverlapAdd(:,i), '-.')
    if i == 1
        semilogy(Fr, CdBtmOverlapAdd(:,i), ':')
    end
    l = legend('NN-IIR', 'UTD-LR', 'UTD', 'BTM', 'Location', 'southeast');
    title(l, 'Diffraction Model')
    xlim([10 100])
    ylim([5e5 1e9])
    xlabel('Update rate (Hz)')
    ylabel('FLOPS')
    title(['Diffraction order: ', num2str(diffOrders(i))])
    colororder(gca, color);
end
fontsize(gcf,scale=textScale)

saveas(gcf, [savePath, filesep, 'SceneDiffractionModelCostsSmall'], 'epsc')
saveas(gcf, [savePath, filesep, 'SceneDiffractionModelCostsSmall'], 'svg')

%% Complexity of calculating UTD

T60 = 0.0375;
rate = 20;
N = ceil(Fs ./ rate + T60 * Fs - 1);

numUtdFreq = 1:N;
CcUtdF = utdOnce + numUtdFreq .* utdPerFreq;

CcNNF = CcNN .* ones(size(numUtdFreq));

figure
plot(numUtdFreq, CcUtdF)
hold on
plot(numUtdFreq, CcNNF)
grid on
xlim([0 512])
l = legend('UTD', 'NN', 'Location', 'northwest')
title(l, 'Diffraction Model')
title({'Cost of calculating UTD and NN'})
xlabel('k (Number of frequencies)')
ylabel('Cost (KFLOPS)')
yticklabels(0:50:300)

saveas(gcf, [savePath, filesep, 'UtdCosts.svg'], 'svg')
