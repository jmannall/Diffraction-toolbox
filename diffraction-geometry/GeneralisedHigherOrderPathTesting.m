close all

%% Data
fs = 96e3;
nfft = 4096;
c = 344;
numEdges = 2;

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

restart = false;
generate = false;
if restart
    i = 1;
end
n = i;
if n == 1
    [meanBTM, meanBTMV2, meanBTMV3, meanBTMV4, meanUTD, meanUTDV2] = deal(zeros(numPaths, 1));
end

dtemplate = struct('rS', [], 'rR', [], 'W', [], 'L', [], 'thetaS', [], 'thetaR', [], 'wedgeIndex', []);
if generate && n == 1
    data = repmat(dtemplate, numPaths, 1);
elseif ~generate
    n = 30;
end
for i = n:numPaths
    if generate
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, data(i)] = GenerateNthOrderPath(numEdges, height);
    else
        wedgeIndex = data(i).wedgeIndex';
        thetaS = data(i).thetaS;
        thetaR = data(i).thetaR;
        radiusS = data(i).rS;
        radiusR = data(i).rR;
        W = data(i).W;
        
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);
    end
    % Generate CAD file from source receiver and top down veiw
    % Calculate geometry parameters (or add as output to above function) for
    % 1st order diffraction sections.
    
    % Expand to include z variation. Requires calculating all the apex points.
    
    createPlot = false;
    %% 2nd order BTM
    
    [ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);
    
    m = 100;
    freq = logspace(log10(20), log10(12e3), m);
    for j = 1:m
        idx(j) = find(fvec < freq(j), 1, "last");
    end
    tfmagDiff2 = tfmag.diff2(idx);
    %% BTM Daisy Chains
    
    % BTM with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [apexTfmag, fvec, apexTfcomplex] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data(i), false, createPlot);
    % BTM with virtual sources and receivers set at rS + W and W + rR from the
    % edge and normalised
    [apexTfmagV2, fvec, apexTfcomplexV2] = SingleBTMApexDaisyChainV2(source, receiver, apex, corners, planeCorners, controlparameters, data(i), false, createPlot);
    % BTM with virtual sources and receivers set at apex points
    [apexTfmagV3, fvec, apexTfcomplexV3] = SingleBTMApexDaisyChainV3(source, receiver, apex, corners, planeCorners, controlparameters, data(i), false, createPlot);
    % BTM with virtual sources and receivers set half way between apex points
    [apexTfmagV4, fvec, apexTfcomplexV4] = SingleBTMApexDaisyChainV4(source, receiver, apex, corners, planeCorners, controlparameters, data(i), false, createPlot);
    
    meanBTM(i) = mean((apexTfmag(idx,end) - tfmagDiff2) .^ 2);
    meanBTMV2(i) = mean((apexTfmagV2(idx,end) - tfmagDiff2) .^ 2);
    meanBTMV3(i) = mean((apexTfmagV3(idx,end) - tfmagDiff2) .^ 2);
    meanBTMV4(i) = mean((apexTfmagV4(idx,end) - tfmagDiff2) .^ 2);
    
    %% UTD
    
    phii = 90;
    
    % UTD with virtual sources and receivers set at apex points and normalised
    % 1 / r
    [utdTfmag, fvec, utdTfcomplex] = UTDSingleApexDaisyChain(data(i), phii, controlparameters, false);
    % UTD with virtual sources and receivers set at apex points
    [utdTfmagV2, fvec, utdTfcomplexV2] = UTDSingleApexDaisyChainV2(data(i), phii, controlparameters, false);
    
    idx = max(2, idx);
    meanUTD(i) = mean((utdTfmag(idx,end) - tfmagDiff2) .^ 2);
    meanUTDV2(i) = mean((utdTfmagV2(idx,end) - tfmagDiff2) .^ 2);
    
    %% Figures
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
%     
%     nexttile([1 3])
%     semilogx(fvec, apexTfmag(:,end))
%     hold on
%     semilogx(fvec, apexTfmagV2(:,end))
%     semilogx(fvec, apexTfmagV3(:,end))
%     semilogx(fvec, apexTfmagV4(:,end))
%     semilogx(fvec, utdTfmag(:,end))
%     semilogx(fvec, utdTfmagV2(:,end))
%     semilogx(fvec, tfmag.diff2)
%     title('Frequency responses')
%     xlabel('Frequency')
%     ylabel('Magnitude')
%     legend('BTM Apex Norm', 'BTM Apex Ext', 'BTM Apex', 'BTM Plane', 'UTD Apex Norm', 'UTD Apex', 'True BTM', 'Location', 'southwest')
%     xlim([20 20000])
%     ylim([-70 0])
end
%% Figures

meanBTMTot = mean(sqrt(meanBTM));
meanBTMV2Tot = mean(sqrt(meanBTMV2));
meanBTMV3Tot = mean(sqrt(meanBTMV3));
meanBTMV4Tot = mean(sqrt(meanBTMV4));
meanUTDTot = mean(sqrt(meanUTD));
meanUTDV2Tot = mean(sqrt(meanUTDV2));

count = 100 * sum(meanBTMV3 > meanBTMV4) / numPaths;
disp('yest')
[N, edges, bin] = histcounts(sqrt(meanBTM), 50);
figure
histogram(sqrt(meanBTM), length(N), 'BinEdges',edges)
title('BTM Apex Norm')
ylim([0 30])

figure
histogram(sqrt(meanBTMV2), length(N), 'BinEdges',edges)
title('BTM Apex Ext')
ylim([0 30])

figure
histogram(sqrt(meanBTMV3), length(N), 'BinEdges',edges)
title('BTM Apex')
ylim([0 30])

figure
histogram(sqrt(meanBTMV4), length(N), 'BinEdges',edges)
title('BTM Plane')
ylim([0 30])

figure
histogram(sqrt(meanUTD), length(N), 'BinEdges',edges)
title('UTD Apex')
ylim([0 30])

figure
histogram(sqrt(meanUTDV2), length(N), 'BinEdges',edges)
title('UTD Apex / 2')
ylim([0 30])
return
%% BTM Daisy Chains

% withCorrection = true;
% [apexTfmag, fvec, apexTfcomplex] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, withCorrection, createPlot);
% [planeApexTfmag, fvec, planeApexTfcomplex] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, false, withCorrection, createPlot);
% [planeTfmag, fvec, planeTfcomplex] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, true, withCorrection, createPlot);

withCorrection = false;
[apexTfmag2, fvec, apexTfcomplex2] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, withCorrection, createPlot);
[apexTfmag2V2, fvec, apexTfcomplex2V2] = SingleBTMApexDaisyChainV2(source, receiver, apex, corners, planeCorners, controlparameters, data, withCorrection, createPlot);

[planeApexTfmag2, fvec, planeApexTfcomplex2] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, false, withCorrection, createPlot);
[planeTfmag2, fvec, planeTfcomplex2] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, true, withCorrection, createPlot);

%% UTD

always = true;
phii = asin(1);
c = 344;

withCorrection = false;
[utdTfmag, fvec, utdTfcomplex] = SingleUTDApexDaisyChain(data, phii, always, c, controlparameters, withCorrection);

% c = 344;
% wedgeIndex = 340;
% v = pi / deg2rad(wedgeIndex);
% thetaS = 10;
% radiusS = 1;
% thetaR = 300;
% radiusR = 1;
% 
% 
% 
% [tfmagUTD, fvec, tfcomplexUTD] = SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phi, always, controlparameters, c);
% 
% [ir, tfmag] = SingleWedge(20, wedgeIndex, thetaS, thetaR, radiusS, radiusR, 10, 10, controlparameters, false);
% 

%% Figures
% figure
% semilogx(fvec, tfmag.diff1)
% hold on
% semilogx(fvec, tfmagUTD)
% legend('BTM', 'UTD')


% figure('Position',[100, 50, 800, 400])
% semilogx(fvec, apexTfmag(:,end))
% hold on
% semilogx(fvec, planeTfmag(:,end))
% semilogx(fvec, planeApexTfmag(:,end))
% semilogx(fvec, utdTfmag(:,end))
% semilogx(fvec, tfmag.diff2)
% title('With correction')
% legend('BTM Apex Daisy Chain', 'BTM Plane Apex Daisy Chain', 'BTM Plane Daisy Chain', 'UTD', 'True BTM', 'Location', 'southwest')
% xlim([20 20000])
% ylim([-70 0])

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
semilogx(fvec, apexTfmag2(:,end))
hold on
semilogx(fvec, apexTfmag2V2(:,end))

semilogx(fvec, planeTfmag2(:,end))
semilogx(fvec, planeApexTfmag2(:,end))
semilogx(fvec, utdTfmag(:,end))
semilogx(fvec, tfmag.diff2)
title('Frequency responses')
xlabel('Frequency')
ylabel('Magnitude')
legend('BTM Apex Daisy Chain', 'BTM Apex Daisy Chain V2', 'BTM Plane Apex Daisy Chain', 'BTM Plane Daisy Chain', 'UTD', 'True BTM', 'Location', 'southwest')
xlim([20 20000])
ylim([-70 0])

%% Test
scale = KimCorrection(data, 2, true);

meanSize = (data.rS + mean(data.W) + data.rR) / 3;
minSize = min(data.rS, min(data.rR, min(data.W)));
maxSize = max(data.rS, max(data.rR, max(data.W)));

