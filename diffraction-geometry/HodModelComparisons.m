close all
% clear all

%% Data
fs = 96e3;
nfft = 8192*2;
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

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges, 'c', c, 'saveFiles', true);

%% Generate geometry

% data = struct('rS', radiusS, 'rR', radiusR, 'W', W, 'L', radiusS + sum(W) + radiusR, 'thetaS', thetaS, 'thetaR', thetaR, 'wedgeIndex', wedgeIndex);
% [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);

numPaths = 100;

index = DataHash({numPaths, fs, nfft, numEdges, height});
[loadPath, savePath] = deal(['geometry/NthOrderPaths_3mTo10m_', num2str(index), '.mat']);
restart = false;
generate = true;
plotFigures = false;
createPlot = false;
if restart
    i = 1;
end
n = i;
%% BTM is very wrong for i = 9. Issue with extra visible paths increasing level and causing comb filtering. - Need to fix and responsible for edge case extreme errors in BTM
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
    [meanBtmApex, meanBtmExt, meanBtmPlane, meanUtdApex, meanUtdExt] = deal(zeros(numPaths, 1));
end

disp('Start')
for i = n:numPaths
    if generate
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, data(i)] = GenerateNthOrderPath(numEdges, height);
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
    [btmTfmagPlane, fvecBTM, btmTfcomplexPlane] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i));
    meanBtmApex(i) = mean((btmTfmagApex(idx,end) - tfmagDiff2) .^ 2);
    meanBtmExt(i) = mean((btmTfmagExt(idx,end) - tfmagDiff2) .^ 2);
    meanBtmPlane(i) = mean((btmTfmagPlane(idx,end) - tfmagDiff2) .^ 2);

    %% UTD
    
    phii = 90;
    
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [utdTfmagApex, ~, utdTfcomplexApex] = SingleUTDApexDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r with Kim correction
    [utdTfmagExt, fvecUTD, utdTfcomplexExt] = SingleUTDExtDaisyChain(data(i), phii, controlparameters);

    x = log10([fvecBTM(2), fvecUTD, fvecBTM(end)]);
    xq = log10(fvecBTM);

    input = [1; zeros(1e3, 1)];
    updateRate = 100;
    windowLength = 2 * fs / updateRate;
    pathLength = data(i).rS + sum(data(i).W) + data(i).rR;
    [~, utdTfmagApexLR] = DelayLineLRFilter(input, utdTfmagApex(:,end), pathLength, windowLength, 1, c, fs, nfft);
    [~, utdTfmagExtLR] = DelayLineLRFilter(input, utdTfmagExt(:,end), pathLength, windowLength, 1, c, fs, nfft);
    for k = 1:numEdges + 1
        utdTfmagApexInterp(:,k) = interp1(x, [utdTfmagApex(1,k); utdTfmagApex(:,k); utdTfmagApex(end,k)], xq);
        utdTfmagExtInterp(:,k) = interp1(x, [utdTfmagExt(1,k); utdTfmagExt(:,k); utdTfmagExt(end,k)], xq);
    end

    xq = log10(fvecBTM(idx));
    xq(1:3) = log10(fvecBTM(2));

    for k = 1:numEdges + 1
        utdTfmagApex2(:,k) = interp1(x, [utdTfmagApex(1,k); utdTfmagApex(:,k); utdTfmagApex(end,k)], xq);
        utdTfmagExt2(:,k) = interp1(x, [utdTfmagExt(1,k); utdTfmagExt(:,k); utdTfmagExt(end,k)], xq);
    end

    meanUtdApex(i) = mean((utdTfmagApexLR(idx) - tfmagDiff2) .^ 2);
    meanUtdExt(i) = mean((utdTfmagExtLR(idx) - tfmagDiff2) .^ 2);
    
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
        semilogx(fvecBTM, btmTfmagApex(:,end))
        hold on
        semilogx(fvecBTM, btmTfmagExt(:,end))
        semilogx(fvecBTM, utdTfmagApexLR)
        semilogx(fvecBTM, utdTfmagExtLR)
        semilogx(fvecBTM, tfmag.diff2)
        title('Frequency responses')
        xlabel('Frequency')
        ylabel('Magnitude')
        legend('BTM Apex', 'BTM Ext', 'UTD Apex', 'UTD Ext', 'True BTM', 'Location', 'southwest')
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
meanUtdApexTot = mean(sqrt(meanUtdApex));
meanUtdExtTot = mean(sqrt(meanUtdExt));

countBtmExtApex = 100 * sum(meanBtmExt < meanBtmApex) / numPaths;
countBtmUtdApex = 100 * sum(meanBtmExt < meanUtdApex) / numPaths;
countBtmUtdExt = 100 * sum(meanBtmExt < meanUtdExt) / numPaths;
countBtmExtPlane = 100 * sum(meanBtmExt < meanBtmPlane) / numPaths;

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
histogram(sqrt(meanBtmPlane), length(N), 'BinEdges',edges)
title('BTM Plane')
ylim([0 30])

figure
histogram(sqrt(meanUtdApex), length(N), 'BinEdges',edges)
title('UTD Apex')
ylim([0 30])

figure
histogram(sqrt(meanUtdExt), length(N), 'BinEdges',edges)
title('UTD Ext')
ylim([0 30])
