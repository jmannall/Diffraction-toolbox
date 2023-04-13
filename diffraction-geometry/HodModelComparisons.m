close all
clear all

set(0, 'DefaultLineLineWidth', 1.5);
colorStore = colororder;

store = [20, 10, 5, 2.5, 1];
for k = 1:length(store)
    %% Data
    fs = 48e3;
    nfft = 8192;
    c = 344;
    numEdges = 2;
    
    height = store(k);
    [zS, zR] = deal(height / 2);
    
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges, 'c', c, 'saveFiles', 3, 'noDirect', true, 'interpolated', false);
    
    %% Generate geometry and controls
    
    numPaths = 500;
    
    index = DataHash({numPaths, fs, nfft, numEdges, 20});
    index = '2180cfff4ddad1de505214480a4535ce';
    [loadPath, savePath] = deal(['geometry/NthOrderPaths_01mTo3m_', num2str(index), '.mat']);
    restart = true;
    generate = false;
    plotFigures = true;
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
    
            wedgeIndex = data(i).wedgeIndex;
            thetaS = data(i).thetaS;
            thetaR = data(i).thetaR;
            radiusS = data(i).rS;
            radiusR = data(i).rR;
            W = data(i).W;
            data(i).L = radiusS + sum(W) + radiusR;
        else
            wI = data(i).wedgeIndex;
            thetaS = data(i).thetaS;
            thetaR = data(i).thetaR;
            rS = data(i).rS;
            rR = data(i).rR;
            W = data(i).W;
            data(i).L = rS + sum(W) + rR;
            %bA = deg2rad(thetaR - thetaS);
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

        X = [deg2rad(wI); bA; mA; wL; [rS, rS + W]; [W + rR, rR]; [zS, zS]; [zR, zR]];
        X = dlarray(single(X), "CB");
    
        [tfmag.NNE(:,1:numEdges), ~, tf.NNE(:,1:numEdges)] = MakeNNPrediction(net, X, pathLength, numFilters, fidx, controlparameters);
    
        tf.NNE(:,end) = 0.5^(numEdges - 1) * (1 / data(i).L) * prod(tf.NNE(:,1:numEdges), 2);
        tfmag.NNE(:,end) = mag2db(abs(tf.NNE(:,end)));
    
        tfmagN.NNE(:,i) = CreateNBandMagnitude(tfmag.NNE(:,end), fidx);

        disp('NN Apex')

        X = [deg2rad(wI); bA; mA; wL; [rS, W]; [W, rR]; [zS, zS]; [zR, zR]];
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
            %tiledlayout(1,2);
    
            %nexttile
            semilogx(fvecBtm, tfmag.Btm, 'LineWidth', 2)
            hold on
            semilogx(fvecBtm, tfmag.BtmE(:,end))
            semilogx(fvec, tfmag.UtdLRE(:,end))
            semilogx(fvec, tfmag.NNE(:,end), '--')
            semilogx(fvecBtm, tfmag.BtmIE(:,end), '--')
            semilogx(fvec, tfmag.UtdILRE(:,end), '--')
            semilogx(fvecBtm, tfmag.BtmIA(:,end), '-.')
            semilogx(fvec, tfmag.UtdLRA(:,end), '-.')
            semilogx(fvec, tfmag.NNA(:,end), ':')
            semilogx(fvecBtm, tfmag.BtmIA(:,end), ':')
            semilogx(fvec, tfmag.UtdILRA(:,end), ':')
            grid on
            colororder(gca, color)
            title('Frequency responses')
            xlabel('Frequency')
            ylabel('Magnitude')
            legend('True BTM', 'BTM Ext', 'UTD Ext', 'NN Ext', 'BTMI Ext', 'UTDI Ext', 'BTM Apex', 'UTD Apex', 'NN Apex', 'BTM ApexI', 'UTD ApexI')
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
        end
    end

%%
    loss = CalculateLoss(tfmagN, tfmagN.Btm);

    %%
    savePath = [savePath 'Test'];
    
    save(savePath, 'data')
    
    disp('Complete')
    
    %% Data
    
    %close all
    
    L = [data.L];
    W = [data.W];
    fields = fieldnames(loss.i);
    numFields = length(fields);

    width = 0.1;
    numBins = 3 / width;
    x = (0:width:3 + width) - width / 2;
    % x = zeros(1, 2 * numBins + 1);
    for i = 2:numBins+1
        idx = width * (i - 2) < W & W <= width * (i - 1);
        for j = 1:numFields
            field = fields{j};
            loss.w.(field)(i) = mean(loss.i.(field)(idx));
        end
    end
    % x = x(1:end - 1);
    for j = 1:numFields
        field = fields{j};
        idx = isnan(loss.w.(field));
        loss.w.(field)(idx) = 0;
        loss.w.(field)(end + 1) = 0;
    end
    
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
    
    save(['hod workspace_Test1deg_', num2str(k)])
end
%%


color = colorStore(1:3,:);
figure
plot(x, loss.w.BtmE)
hold on
grid on
plot(x, loss.w.UtdLRE)
plot(x, loss.w.NNE)
plot(x, loss.w.BtmA, '--')
plot(x, loss.w.UtdLRA, '--')
plot(x, loss.w.NNA, '--')
plot(x, loss.w.BtmIE, '-.')
plot(x, loss.w.UtdILRE, '-.')
plot(x, loss.w.NNE)
plot(x, loss.w.BtmIA, ':')
plot(x, loss.w.UtdILRA, ':')
colororder(color)
xlim([0.1 3])
legend('BTM Extension', 'UTD-LR Extension', 'NN-IIR Extension', 'BTM Apex', 'UTD-LR Apex', 'NN-IIR Apex', 'BTMI Extension', 'UTDI-LR Extension', '-', 'BTMI Apex', 'UTDI-LR Apex')
xlabel('W_1 (m)')
ylabel('Mean absolute error (dBA)')

saveDir = 'figures';
saveas(gcf, [saveDir filesep 'HODComparison_Test'], 'epsc')
saveas(gcf, [saveDir filesep 'HODComparison_Test'], 'svg')

%% Result Plot

%close all

%load(['hod workspace_Test1deg_' num2str(k), '.mat'])
color = colorStore([4,1,2,5,6,3],:);

figure
plot(x, loss.w.BtmIE)
hold on
grid on
plot(x, loss.w.NNE)
plot(x, loss.w.UtdILRE)
plot(x, loss.w.BtmIA, '--')
plot(x, loss.w.NNA, '--')
plot(x, loss.w.UtdILRA, '--')
%plot(x, meanUEI)
colororder(color)
xlim([0.1 3])
ylim([0 6])
title(num2str(store(k)))
legend('BTM Extension', 'NN-IIR Extension', 'UTD-LR Extension', 'BTM Apex', 'NN-IIR Apex', 'UTD-LR Apex')
xlabel('W_1 (m)')
ylabel('Mean absolute error (dBA)')


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