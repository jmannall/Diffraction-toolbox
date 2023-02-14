%close all
%clear all
set(0, 'DefaultLineLineWidth', 1.5);
colorStore = colororder;

store = [20, 10, 5, 2.5, 1];
for k = 1:length(store)
%% Data
fs = 96e3;
nfft = 8192*2;
c = 344;
numEdges = 2;
fvec = fs / nfft * (0:nfft / 2 - 1);

% thetaS = 10;
% thetaR = 205;
% wedgeIndex = [220; 240];
% W = 1.6;
% radiusS = 1.1;
% radiusR = 0.9;
height = store(k);
[zS, zR] = deal(height / 2);

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges, 'c', c, 'saveFiles', 3, 'noDirect', true, 'interpolated', false);

%% Generate geometry

% data = struct('rS', radiusS, 'rR', radiusR, 'W', W, 'L', radiusS + sum(W) + radiusR, 'thetaS', thetaS, 'thetaR', thetaR, 'wedgeIndex', wedgeIndex);
% [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);

numPaths = 500;

index = DataHash({numPaths, fs, nfft, numEdges, 20});
[loadPath, savePath] = deal(['geometry/NthOrderPaths_01mTo3m_', num2str(index), '.mat']);
restart = true;
generate = false;
plotFigures = false;
createPlot = false;
if restart
    i = 1;
end
n = i;
%% BTM is very wrong for i = 9. Issue with extra visible paths increasing level and causing comb filtering. - Need to fix and responsible for edge case extreme errors in BTM

dtemplate = struct('rS', [], 'rR', [], 'W', [], 'L', [], 'thetaS', [], 'thetaR', [], 'wedgeIndex', []);
if generate && n == 1
    data = repmat(dtemplate, numPaths, 1);
elseif ~generate
    %n = 1;
    load(loadPath, 'data')
end

if n == 1
    %[meanBtmApex, meanBtmExt, meanBtmPlane, meanNNExt, meanNNApex, meanUtdApex, meanUtdExt, meanBtmApexI, meanBtmExtI, meanUtdApexI, meanUtdExtI, meanUtdExtI1, meanNNExt1] = deal(zeros(numPaths, 1));
end

[~, tfmagDefault, ~, fvec] = DefaultBTM(controlparameters);
octBand = 8;
tfmagDefault = tfmagDefault(1:end / 2,:);

fs = controlparameters.fs / 2;
nfft = controlparameters.nfft / 2;
fvec = fs/nfft*[0:nfft/2-1];
[targetData, fc, fidx] = CreateFrequencyNBands(tfmagDefault, fvec, octBand);

AWeighting = weightingFilter('A-weighting',fs);
aWeight = freqz(AWeighting, fc, fs);
aWeight = abs(aWeight)';

% for i = 1:numPaths
%     data(i).rS = 10;
%     data(i).rR = 15;
%     data(i).W = i / 20;
%     data(i).thetaS = 30;
%     data(i).thetaR = 225;
% end

%% NN Prep
    
loadDir = 'NNSaves';

netName = '5_iir-2057_0001-1-09-099-3-25.mat';

numFilters = 2;
filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
load([loadDir, filesep, netName], 'net');

wedgeLength = [height, height];
biquad = false;

windowLength = 12;
input = [1; zeros(23, 1)];

disp('Start')
for i = n:numPaths
    if generate
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, data(i)] = GenerateNthOrderPath(numEdges, height);

        wedgeIndex = data(i).wedgeIndex;
        thetaS = data(i).thetaS;
        thetaR = data(i).thetaR;
        radiusS = data(i).rS;
        radiusR = data(i).rR;
        W = data(i).W;
        data(i).L = radiusS + sum(W) + radiusR;
    else
%         height = 3;
%         data(i).wedgeIndex = [200 200 200 200 200];
%         data(i).thetaS = 5;
%         data(i).thetaR = 195;
%         data(i).rS = 0.5;
%         data(i).rR = 0.5;
%         data(i).W = [0.2 0.2 0.2 0.2];
%         controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 5, 'c', c);

%         height = 3;
%         data(i).wedgeIndex = [270 270];
%         data(i).thetaS = 85;
%         data(i).thetaR = 185;
%         data(i).rS = 1.4;
%         data(i).rR = 1.4;
%         data(i).W = 0.5 / 3.5;

        wedgeIndex = data(i).wedgeIndex;
        thetaS = data(i).thetaS;
        thetaR = data(i).thetaR;
        radiusS = data(i).rS;
        radiusR = data(i).rR;
        W = data(i).W;
        data(i).L = radiusS + sum(W) + radiusR;

        [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid, vReceiver] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);
    end

%     zS = 0.1 + 4.8 * rand(1);   % Not quite working as BtmDaisy chain nt set up to work right with z values
%     zR = 0.1 + 4.8 * rand(1);
% 
%     data(i).zS = zS;
%     data(i).zR = zR;
% 
%     source(:,3) = zS;    
%     receiver(:,3) = zR;

    % Expand to include z variation. Requires calculating all the apex points.
    
    %% 2nd order BTM
    %PlotGeometry(corners, planeCorners, source, receiver)

    [ir, tfmag, tvec, fvecBTM, tfcomplex] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);
    
%     vSource = [source; apex(1:numEdges - 1,:)];
%     vReceiver = [apex(2:numEdges,:); receiver];
%     rS = [data(i).rS, data(i).W];
%     rR = [data(i).rR, fliplr(data(i).W)];
% 
%     cumRs = cumsum(rS)';
%     cumRr = fliplr(cumsum(rR))';
% 
%     vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
%     vector = vSource - vCorners;
%     vector = cumRs .* vector ./ vecnorm(vector, 2, 2);
%     vSource = vCorners + vector;
% 
%     vector = vReceiver - vCorners;
%     vector = cumRr .* vector ./ vecnorm(vector, 2, 2);
%     vReceiver = vCorners + vector;
% 
%     epsilon = 1e-2;

%     test = [apex; vReceiver{1}];
%     test2 = [apex; vReceiver{2}];
%     figure
%     plot3(source(:,1), source(:,2), source(:,3), 'o')
%     hold on
%     plot3(receiver(:,1), receiver(:,2), receiver(:,3), 'o')
%     plot3(test(:,1), test(:,2), test(:,3))
%     plot3(test2(:,1), test2(:,2), test2(:,3), '--')
%     plot3(Q(:,1), Q(:,2), Q(:,3))
%     plot3(apex(:,1), apex(:,2), apex(:,3), 'o')
%     title('Scene')
%     legend('Source', 'Receiver', 'vR', 'vR', 'Planes', 'Apex')
%     view([0 90])
%     xlim([-4 2])
%     ylim([-4 2])
% 
%     [~, tfmagDiff1Ref, ~, ~, tfcomplexDiff1Ref] = SingleBTM(source, vReceiver{1}, corners, planeCorners, planeRigid, controlparameters, createPlot);
%     [~, tfmagDiff2Ref, ~, ~, tfcomplexDiff2Ref] = SingleBTM(source, vReceiver{2}, corners, planeCorners, planeRigid, controlparameters, createPlot);

    tfmagDiff2 = CreateNBandMagnitude(tfmag.diff2, fidx);
    %% BTM Daisy Chains
    
    % BTM with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [btmTfmagApex, ~, btmTfcomplexApex] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    % BTM with virtual sources and receivers set at rS + W and W + rR from the
    % edge and normalised
    [btmTfmagExt, ~, btmTfcomplexExt] = SingleBTMExtDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    
    btmTfmagApexNBand = CreateNBandMagnitude(btmTfmagApex(:,end), fidx);
    btmTfmagExtNBand = CreateNBandMagnitude(btmTfmagExt(:,end), fidx);

    %% BTM Daisy Chains interpolated

    epsilon = 1e-2;
    [btmTfmagExtI, btmTfcomplexExtI, btmTfmagApexI, btmTfcomplexApexI] = deal(zeros(nfft, numEdges + 1));
    [btmTfmagExtI(:,1), ~, btmTfcomplexExtI(:,1)] = SingleWedgeInterpolated(wedgeLength(1), wedgeIndex(1), epsilon, wedgeIndex(1) - thetaS, radiusS, W + radiusR, zS, zR, controlparameters, false);
    [btmTfmagExtI(:,2), ~, btmTfcomplexExtI(:,2)] = SingleWedgeInterpolated(wedgeLength(2), wedgeIndex(2), epsilon, thetaR, radiusS + W, radiusR, zS, zR, controlparameters, false);
    [btmTfmagApexI(:,1), ~, btmTfcomplexApexI(:,1)] = SingleWedgeInterpolated(wedgeLength(1), wedgeIndex(1), epsilon, wedgeIndex(1) - thetaS, radiusS, W, zS, zR, controlparameters, false);
    [btmTfmagApexI(:,2), ~, btmTfcomplexApexI(:,2)] = SingleWedgeInterpolated(wedgeLength(2), wedgeIndex(2), epsilon, thetaR, W, radiusR, zS, zR, controlparameters, false);

    btmTfcomplexExtI(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(btmTfcomplexExtI(:,1:numEdges), 2);
    btmTfcomplexApexI(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(btmTfcomplexApexI(:,1:numEdges), 2);
    btmTfmagExtI(:,end) = mag2db(abs(btmTfcomplexExtI(:,end)));
    btmTfmagApexI(:,end) = mag2db(abs(btmTfcomplexApexI(:,end)));

    btmTfmagExtINBand = CreateNBandMagnitude(btmTfmagExtI(:,end), fidx);
    btmTfmagApexINBand = CreateNBandMagnitude(btmTfmagApexI(:,end), fidx);

    %% NN Ext 

    disp('NN Ext')
    output = CreateNNOutput(net, wedgeIndex', wedgeLength', [wedgeIndex(1) - epsilon; thetaR], [thetaS; epsilon], [radiusS; radiusS + W], [W + radiusR; radiusR], [zS; zS], [zR; zR], false);
    pathLength = [data(i).L, data(i).L];
    validPath = true(1, numEdges);
    [tfcomplexNNExt, b, a] = CalculateNN(output, btmTfcomplexExt(:,1:numEdges), validPath, [1;1], nfft, fs, biquad);

    [bAll, aAll] = deal(zeros(2, 2 * numEdges));
    for j = 1:numEdges
        idx = 2 * j - 1:2 * j;
        bAll(:, idx) = b(:,:,j);
        aAll(:, idx) = a(:,:,j);
    end

    [numCoeff, numFilters, ~] = size(bAll);
    [tempInputBufferIIR, tempOutputBufferIIR] = InitialiseIIRBuffers(numCoeff, numFilters);
    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs, true);
    amplitude = 0.5 .^ (numEdges - 1) * amplitude;
    iirLength = 1e4;
    irNN = zeros(delay + iirLength + 1, 1);
    irIn = zeros(delay + iirLength + 1, 1);
    irIn(delay) = amplitude * (1 - fracDelay);
    irIn(delay + 1) = amplitude * fracDelay;
    for j = delay:delay + iirLength
        [tempInputBufferIIR, tempOutputBufferIIR, irNN(j)] = ProcessIIRFilter(irIn(j), bAll, aAll, tempInputBufferIIR, tempOutputBufferIIR);
    end
    [~, tfcomplexNNExt(:,3)] = IrToTf(irNN, nfft);

    tfmagNNExt = mag2db(abs(tfcomplexNNExt));
    tfcomplexNN = 0.5^(numEdges - 1) * 1 / pathLength(1) .* prod(tfcomplexNNExt(:,1:numEdges),2);

    tfmagNN = mag2db(abs(tfcomplexNN));

    tfmagNNExtNBand = CreateNBandMagnitude(tfmagNNExt(:,end), fidx);

    %% NN Apex

    disp('NN Apex')
    output = CreateNNOutput(net, wedgeIndex', wedgeLength', [wedgeIndex(1) - epsilon; thetaR], [thetaS; epsilon], [radiusS; W], [W; radiusR], [zS; zS], [zR; zR], false);
    pathLength = [data(i).L, data(i).L];
    validPath = true(1, numEdges);
    [tfcomplexNNApex, b, a] = CalculateNN(output, btmTfcomplexExt(:,1:numEdges), validPath, [1;1], nfft, fs, biquad);

    [bAll, aAll] = deal(zeros(2, 2 * numEdges));
    for j = 1:numEdges
        idx = 2 * j - 1:2 * j;
        bAll(:, idx) = b(:,:,j);
        aAll(:, idx) = a(:,:,j);
    end

    [numCoeff, numFilters, ~] = size(bAll);
    [tempInputBufferIIR, tempOutputBufferIIR] = InitialiseIIRBuffers(numCoeff, numFilters);
    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs, true);
    amplitude = 0.5 .^ (numEdges - 1) * amplitude;
    iirLength = 1e4;
    irNN = zeros(delay + iirLength + 1, 1);
    irIn = zeros(delay + iirLength + 1, 1);
    irIn(delay) = amplitude * (1 - fracDelay);
    irIn(delay + 1) = amplitude * fracDelay;
    for j = delay:delay + iirLength
        [tempInputBufferIIR, tempOutputBufferIIR, irNN(j)] = ProcessIIRFilter(irIn(j), bAll, aAll, tempInputBufferIIR, tempOutputBufferIIR);
    end
    [~, tfcomplexNNApex(:,3)] = IrToTf(irNN, nfft);
    
    tfmagNNApex = mag2db(abs(tfcomplexNNApex(:,3)));

    tfmagNNApexNBand = CreateNBandMagnitude(tfmagNNApex(:,end), fidx);

%     output = CreateNNOutput(net, wedgeIndex', wedgeLength', [wedgeIndex(1) - 0.001; thetaR], [thetaS; 0.001], [radiusS; W], [W; radiusR], [zS; zS], [zR; zR], false);
%     pathLength = [data(i).L, data(i).L];
%     validPath = true(1, numEdges);
%     [tfcomplexNNApex, b, a] = CalculateNN(output, btmTfcomplexExt(:,1:numEdges), validPath, [1;1], nfft, fs, biquad);
% 
%     [~, irNN] = DelayLineIIRFilter(input, pathLength, windowLength, validPath, b, a, c, fs / 2, true);
%     [~, tfcomplexNNApex(:,3)] = IrToTf(irNN, nfft);
% 
%     tfcomplexNNApex(:,3) = 0.5^(numEdges - 1) * 1 / pathLength(1) .* prod(tfcomplexNNApex(:,1:numEdges),2);
%     tfmagNNApex = mag2db(abs(tfcomplexNNApex));
% 
%     tfmagNNApexNBand = CreateNBandMagnitude(tfmagNNApex(:,end), fidx);
%     meanNNApex(i) = mean(aWeight .* abs(tfmagNNApexNBand - tfmagDiff2));

    %% UTD

    controlparameters.interpolated = false;
    
    phii = 90;
    %[zA, phii] = CalculateApex(W + radiusR, radiusS, zS, zR, wedgeLength, true);
    
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [utdTfmagApex, ~, utdTfcomplexApexI] = SingleUTDApexDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r with Kim correction

    %[zA, phii] = CalculateApex([W + radiusR, radiusR], [radiusS, radiusS + W], [zS, zS], [zR, zR], wedgeLength, true);
    [utdTfmagExt, fvecUTD, utdTfcomplexExtI] = SingleUTDExtDaisyChain(data(i), phii, controlparameters);

    input = [1; zeros(1e3, 1)];
    updateRate = 100;
    windowLength = 2 * fs / updateRate;
    pathLength = data(i).rS + sum(data(i).W) + data(i).rR;
    [utdTfmagApexLR, utdTfmagExtLR]    = deal(zeros(nfft / 2, numEdges));
    for j = 1:numEdges + 1
        [~, utdIrApexLR] = DelayLineLR(input, pathLength, windowLength, 1, utdTfmagApex(:,j)', c, fs, 1);
        [~, utdIrExtLR] = DelayLineLR(input, pathLength, windowLength, 1, utdTfmagExt(:,j)', c, fs, 1);
        utdTfmagApexLR(:,j) = IrToTf(utdIrApexLR, nfft);
        utdTfmagExtLR(:,j) = IrToTf(utdIrExtLR, nfft);
    end

    utdTfmagApexLRNBand = CreateNBandMagnitude(utdTfmagApexLR(:,end), fidx);
    utdTfmagExtLRNBand = CreateNBandMagnitude(utdTfmagExtLR(:,end), fidx);
    
    %% UTD Interpolated

    controlparameters.interpolated = true;

    %phii = 90;
    [zA, phii] = CalculateApex(W + radiusR, radiusS, zS, zR, wedgeLength, true);
    
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [utdTfmagApexI, ~, utdTfcomplexApexI] = SingleUTDApexDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r with Kim correction

    %[zA, phii] = CalculateApex([W + radiusR, radiusR], [radiusS, radiusS + W], [zS, zS], [zR, zR], wedgeLength, true);
    [utdTfmagExtI, ~, utdTfcomplexExtI] = SingleUTDExtDaisyChain(data(i), phii(1), controlparameters);

    input = [1; zeros(1e3, 1)];
    updateRate = 100;
    windowLength = 2 * fs / updateRate;
    pathLength = data(i).rS + sum(data(i).W) + data(i).rR;
    [utdTfmagApexILR, utdTfmagExtILR]    = deal(zeros(nfft / 2, numEdges));
    for j = 1:numEdges + 1
        [~, utdIIrApexLR] = DelayLineLR(input, pathLength, windowLength, 1, utdTfmagApexI(:,j)', c, fs, 1);
        [~, utdIIrExtLR] = DelayLineLR(input, pathLength, windowLength, 1, utdTfmagExtI(:,j)', c, fs, 1);
        utdTfmagApexILR(:,j) = IrToTf(utdIIrApexLR, nfft);
        utdTfmagExtILR(:,j) = IrToTf(utdIIrExtLR, nfft);
    end

    utdTfmagApexILRNBand = CreateNBandMagnitude(utdTfmagApexILR(:,end), fidx);
    utdTfmagExtILRNBand = CreateNBandMagnitude(utdTfmagExtILR(:,end), fidx);

    %% 1st order N Bands

    btmTfmagExtINBand1 = CreateNBandMagnitude(btmTfmagExtI(:,1:2), fidx);
    utdTfmagExtILRNBand1 = CreateNBandMagnitude(utdTfmagExtILR(:,1:2), fidx);
    tfmagNNExtNBand1 = CreateNBandMagnitude(tfmagNNExt(:,1:2), fidx);

    %% Other mean figures

    meanBtmExt(i) = mean(aWeight .* abs(btmTfmagExtNBand - tfmagDiff2));
    meanUtdExt(i) = mean(aWeight .* abs(utdTfmagExtLRNBand - tfmagDiff2));

    meanBtmApex(i) = mean(aWeight .* abs(btmTfmagApexNBand - tfmagDiff2));
    meanUtdApex(i) = mean(aWeight .* abs(utdTfmagApexLRNBand - tfmagDiff2));

    meanBtmExtI(i) = mean(aWeight .* abs(btmTfmagExtINBand - tfmagDiff2));
    meanUtdExtI(i) = mean(aWeight .* abs(utdTfmagExtILRNBand - tfmagDiff2));
    meanNNExt(i) = mean(aWeight .* abs(tfmagNNExtNBand - tfmagDiff2));

    meanBtmApexI(i) = mean(aWeight .* abs(btmTfmagApexINBand - tfmagDiff2));
    meanUtdApexI(i) = mean(aWeight .* abs(utdTfmagApexILRNBand - tfmagDiff2));
    meanNNApex(i) = mean(aWeight .* abs(tfmagNNApexNBand - tfmagDiff2));

    meanUtdExtI1(i) = mean(mean(aWeight .* abs(utdTfmagExtILRNBand1 - btmTfmagExtINBand1)));
    meanNNExt1(i) = mean(mean(aWeight .* abs(tfmagNNExtNBand1 - btmTfmagExtINBand1)));

    meanUtdExtI2(i) = mean(mean(aWeight .* abs(utdTfmagExtILRNBand - btmTfmagExtINBand)));
    meanNNExt2(i) = mean(mean(aWeight .* abs(tfmagNNExtNBand - btmTfmagExtINBand)));

    %% Figures
    if plotFigures
        %close all

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
        figure('Position',[50 100 1820 800])
        tiledlayout(1,2);

        nexttile
        semilogx(fvecBTM, tfmag.diff2, 'LineWidth', 2)
        hold on
        semilogx(fvecBTM, btmTfmagExt(:,end))
        semilogx(fvec, utdTfmagExtLR(:,end))
        semilogx(fvec, tfmagNNExt(:,end), '--')
        semilogx(fvecBTM, btmTfmagExtI(:,end), '--')
        semilogx(fvec, utdTfmagExtILR(:,end), '--')
        semilogx(fvecBTM, btmTfmagApex(:,end), '-.')
        semilogx(fvec, utdTfmagApexLR(:,end), '-.')
        semilogx(fvec, tfmagNNApex(:,end), ':')
        semilogx(fvecBTM, btmTfmagApexI(:,end), ':')
        semilogx(fvec, utdTfmagApexILR(:,end), ':')
        grid on
        colororder(gca, color)
        title('Frequency responses')
        xlabel('Frequency')
        ylabel('Magnitude')
        legend('True BTM', 'BTM Ext', 'UTD Ext', 'NN Ext', 'BTM ExtI', 'UTD ExtI', 'BTM Apex', 'UTD Apex', 'NN Apex', 'BTM ApexI', 'UTD ApexI')
        xlim([20 20000])
        ylim([-70 0])

        color = colorStore(1:5,:);
        nexttile
        for j = 1:2
            semilogx(fvec, tfmagNNExt(:,j), 'LineWidth', 2)
            hold on
            semilogx(fvecBTM, btmTfmagExt(:,j))
            semilogx(fvecBTM, btmTfmagExtI(:,j), '--')
            semilogx(fvec, utdTfmagExtLR(:,j), '-.')
            semilogx(fvec, utdTfmagExtILR(:,j), ':')
        end
        grid on
        colororder(gca, color)
        title('Frequency responses')
        xlabel('Frequency')
        ylabel('Magnitude')
        legend('NN', 'BTM Ext', 'BTM ExtI', 'UTD Ext', 'UTD ExtI')
        xlim([20 20000])
        ylim([-30 10])
    end
end

save(savePath, 'data')

disp('Complete')

%% Data

%close all

L = [data.L];
W = [data.W];

width = 0.1;
numBins = 3 / width;
[meanBA, meanBE, meanNE, meanNA, meanUE, meanUA, meanBAI, meanBEI, meanUAI, meanUEI, meanNN1, meanUTD1, meanNN2, meanUTD2] = deal(zeros(1, numBins + 2));
x = (0:width:3 + width) - width / 2;
% x = zeros(1, 2 * numBins + 1);
for i = 2:numBins+1
    idx = width * (i - 2) < W & W <= width * (i - 1);

    meanBE(i) = mean(meanBtmExt(idx));
    meanUE(i) = mean(meanUtdExt(idx));

    meanBA(i) = mean(meanBtmApex(idx));
    meanUA(i) = mean(meanUtdApex(idx));

    meanBEI(i) = mean(meanBtmExtI(idx));
    meanUEI(i) = mean(meanUtdExtI(idx));
    meanNE(i) = mean(meanNNExt(idx));

    meanBAI(i) = mean(meanBtmApexI(idx));
    meanUAI(i) = mean(meanUtdApexI(idx));
    meanNA(i) = mean(meanNNApex(idx));

    meanNN1(i) = mean(meanNNExt1(idx));
    meanUTD1(i) = mean(meanUtdExtI1(idx));

    meanNN2(i) = mean(meanNNExt2(idx));
    meanUTD2(i) = mean(meanUtdExtI2(idx));
end
% x = x(1:end - 1);
idx = isnan(meanBA);
[meanBE(idx), meanUE(idx)] = deal(0);
[meanBA(idx), meanUA(idx)] = deal(0);
[meanNE(idx), meanBEI(idx), meanUEI(idx)] = deal(0);
[meanNA(idx), meanBAI(idx), meanUAI(idx)] = deal(0);
[meanNN1(idx), meanUTD1(idx), meanNN2(idx), meanUTD2(idx)] = deal(0);

figure
plot(x, meanNN1)
hold on
grid on
plot(x, meanUTD1, '--')
plot(x, meanNN2)
plot(x, meanUTD2, '--')
plot(x, meanNE)
plot(x, meanUEI, '--')
plot(x, meanBEI, '-.')
xlim([0.1 3])
legend('meanNNExt1', 'meanUtdExtI1', 'meanNNExt2', 'meanUtdExtI2', 'meanNNExt', 'meanUtdExtI', 'meanBtmExtI')
title(num2str(height))
xlabel('W_1 (m)')
ylabel('Mean absolute error (dBA)')

saveDir = 'figures';
saveas(gcf, [saveDir filesep 'HODComparison_TestNew_' num2str(k)], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparison_TestNew_' num2str(k)], 'svg')

save(['hod workspace_Test_', num2str(k)])
end


color = colorStore(1:3,:);
figure
plot(x, meanBE)
hold on
grid on
plot(x, meanUE)
plot(x, meanNE)
plot(x, meanBA, '--')
plot(x, meanUA, '--')
plot(x, meanNA, '--')
plot(x, meanBEI, '-.')
plot(x, meanUEI, '-.')
plot(x, meanNE)
plot(x, meanBAI, ':')
plot(x, meanUAI, ':')
colororder(color)
xlim([0.1 3])
legend('BTM Extension', 'UTD-LR Extension', 'NN-IIR Extension', 'BTM Apex', 'UTD-LR Apex', 'NN-IIR Apex', 'BTMI Extension', 'UTDI-LR Extension', '-', 'BTMI Apex', 'UTDI-LR Apex')
xlabel('W_1 (m)')
ylabel('Mean absolute error (dBA)')

saveDir = 'figures';
saveas(gcf, [saveDir filesep 'HODComparison_Test'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparison_Test'], 'svg')

%% Result Plot

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