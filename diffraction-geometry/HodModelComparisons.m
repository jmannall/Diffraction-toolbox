close all
clear all

%% Data
fs = 96e3;
nfft = 4096;
c = 344;
numEdges = 2;
fvec = fs / nfft * (0:nfft / 2 - 1);

thetaS = 10;
thetaR = 205;
wedgeIndex = [220; 240];
W = 1.6;
radiusS = 1.1;
radiusR = 0.9;
height = 20;
[zS, zR] = deal(height / 2);

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges, 'c', c);

%% Generate geometry

% data = struct('rS', radiusS, 'rR', radiusR, 'W', W, 'L', radiusS + sum(W) + radiusR, 'thetaS', thetaS, 'thetaR', thetaR, 'wedgeIndex', wedgeIndex);
% [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);

numPaths = 100;

index = DataHash({numPaths, fs, nfft, numEdges, height});
[loadPath, savePath] = deal(['geometry/NthOrderPaths_', num2str(index), '.mat']);
restart = false;
generate = false;
plotFigures = false;
createPlot = false;
if restart
    i = 1;
end
%n = i;

m = 100;
freq = logspace(log10(20), log10(12e3), m);
idx = zeros(1, m);
for j = 1:m
    idx(j) = find(fvec < freq(j), 1, "last");
end
idxUtd = max(2, idx);

dtemplate = struct('rS', [], 'rR', [], 'W', [], 'L', [], 'thetaS', [], 'thetaR', [], 'wedgeIndex', []);
if generate && n == 1
    data = repmat(dtemplate, numPaths, 1);
elseif ~generate
    n = 1;
    load(loadPath, 'data')
end

if n == 1
    [meanBtmApex, meanBtmExt, meanBtmPlane, meanUtd, meanUtdKim] = deal(zeros(numPaths, 1));
end

disp('Start')
for i = n:numPaths
    if generate
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, data(i)] = GenerateNthOrderPath(numEdges, height);
    else
        wedgeIndex = data(i).wedgeIndex;
        thetaS = data(i).thetaS;
        thetaR = data(i).thetaR;
        radiusS = data(i).rS;
        radiusR = data(i).rR;
        W = data(i).W;
        
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);
    end
    
    % Expand to include z variation. Requires calculating all the apex points.
    
    %% 2nd order BTM
    
    [ir, tfmag, tvec, ~, tfcomplex] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);
    
    tfmagDiff2 = tfmag.diff2(idx);
    %% BTM Daisy Chains
    
    % BTM with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [btmTfmagApex, ~, btmTfcomplexApex] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    % BTM with virtual sources and receivers set at rS + W and W + rR from the
    % edge and normalised
    [btmTfmagExt, ~, btmTfcomplexExt] = SingleBTMExtDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    [btmTfmagPlane, ~, btmTfcomplexPlane] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    meanBtmApex(i) = mean((btmTfmagApex(idx,end) - tfmagDiff2) .^ 2);
    meanBtmExt(i) = mean((btmTfmagExt(idx,end) - tfmagDiff2) .^ 2);
    meanBtmPlane(i) = mean((btmTfmagPlane(idx,end) - tfmagDiff2) .^ 2);

    %% UTD
    
    phii = 90;
    
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [utdTfmag, ~, utdTfcomplex] = UTDSingleDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r with Kim correction
    [utdTfmagKim, fvec, utdTfcomplexKim] = UTDSingleDaisyChain(data(i), phii, controlparameters, true);
    
    meanUtd(i) = mean((utdTfmag(idxUtd,end) - tfmagDiff2) .^ 2);
    meanUtdKim(i) = mean((utdTfmagKim(idxUtd,end) - tfmagDiff2) .^ 2);
    
    %% Figures
    if plotFigures
        figure('Position',[50 100 1820 800])
        t = tiledlayout(1,5);
        nexttile([1 2])
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
        
        nexttile([1 3])
        semilogx(fvec, btmTfmagApex(:,end))
        hold on
        semilogx(fvec, btmTfmagExt(:,end))
        semilogx(fvec, btmTfmagPlane(:,end))
        semilogx(fvec, utdTfmag(:,end))
        semilogx(fvec, utdTfmagKim(:,end))
        semilogx(fvec, tfmag.diff2)
        title('Frequency responses')
        xlabel('Frequency')
        ylabel('Magnitude')
        legend('BTM Apex', 'BTM Ext', 'BTM Plane', 'UTD', 'UTD Kim', 'True BTM', 'Location', 'southwest')
        xlim([20 20000])
        ylim([-70 0])
    end
end

save(savePath, 'data')

disp('Complete')

%% Figures
close all

meanBtmApexTot = mean(sqrt(meanBtmApex));
meanBtmExtTot = mean(sqrt(meanBtmExt));
meanBtmPlaneTot = mean(sqrt(meanBtmPlane));
meanUtdTot = mean(sqrt(meanUtd));
meanUtdKimTot = mean(sqrt(meanUtdKim));

countBtmExtApex = 100 * sum(meanBtmExt < meanBtmApex) / numPaths;
countBtmUtd = 100 * sum(meanBtmExt < meanUtd) / numPaths;
countBtmUtdKim = 100 * sum(meanBtmExt < meanUtdKim) / numPaths;
countBtmExtPlane = 100 * sum(meanBtmExt < meanBtmPlane) / numPaths;

[N, edges, bin] = histcounts(sqrt(meanBtmApex), 50);

figure
histogram(sqrt(meanBtmApex), length(N), 'BinEdges',edges)
title('BTM Apex')
ylim([0 30])

figure
histogram(sqrt(meanBtmExt), length(N), 'BinEdges',edges)
title('BTM Ext')
ylim([0 30])

figure
histogram(sqrt(meanBtmPlane), length(N), 'BinEdges',edges)
title('BTM Plane')
ylim([0 30])

figure
histogram(sqrt(meanUtd), length(N), 'BinEdges',edges)
title('UTD')
ylim([0 30])

figure
histogram(sqrt(meanUtdKim), length(N), 'BinEdges',edges)
title('UTD Kim')
ylim([0 30])
