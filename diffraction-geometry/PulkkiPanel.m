close all
clear all

width = 0.47;
depth = 0.03;
height = 0.62;

fs = 48e3;
nfft = 8192;
c= 344;

controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'saveFiles', 3, 'difforder', 1, 'noDirect', false);
createPlot = true;

[ir, tfmag, tvec, fvec, tfcomplex] = SinglePanel(width, depth, height, controlparameters, createPlot);

thetaRGlobal = (0.1:1:180.1)';
PlotSpectrogram([tfcomplex.direct], fvec, thetaRGlobal', [-15 0], 'Pulkki Panel Direct', false, false, 'Receiver azimuth')
PlotSpectrogram([tfcomplex.geom], fvec, thetaRGlobal', [-15 0], 'Pulkki Panel Specular', false, false, 'Receiver azimuth')
PlotSpectrogram([tfcomplex.diff1], fvec, thetaRGlobal', [-15 0], 'Pulkki Panel Diffraction', false, false, 'Receiver azimuth')
PlotSpectrogram([tfcomplex.complete], fvec, thetaRGlobal', [-15 0], 'Pulkki Panel', false, false, 'Receiver azimuth')

load('SinglePanel_91c6282718b3d3e0fb9222e6ae5ab5b3_Sdata.mat')
load('SinglePanel_91c6282718b3d3e0fb9222e6ae5ab5b3_Rdata.mat')
load('SinglePanel_91c6282718b3d3e0fb9222e6ae5ab5b3_paths.mat')
load('SinglePanel_91c6282718b3d3e0fb9222e6ae5ab5b3_eddata.mat')
%load('SinglePanel_91c6282718b3d3e0fb9222e6ae5ab5b3_ed2data.mat')

%% NN

numEdges = 12;
numReceivers = length(thetaRGlobal);
const = ones(numEdges, numReceivers);

wedgeLength = edgedata.edgelengthvec;
wedgeIndex = 360 - rad2deg(edgedata.closwedangvec);
rS = Sdata.rSsho(Sdata.reftoshortlistS);
zS = Sdata.zSsho(Sdata.reftoshortlistS);
thetaS = rad2deg(Sdata.thetaSsho(Sdata.reftoshortlistS));
rR = Rdata.rRsho(Rdata.reftoshortlistR);
zR = Rdata.zRsho(Rdata.reftoshortlistR);
thetaR = rad2deg(Rdata.thetaRsho(Rdata.reftoshortlistR));

loadDir = 'NNSaves';
netName = 'iir-2057_0001-1-09-099-3-25.mat';
%netName = '5_iir-2057_0001-1-09-099-3-25.mat';

[~, tfmagDefault, ~, fvec] = DefaultBTM(controlparameters);
[tfmagNBand, fc, fidx] = CreateFrequencyNBands([tfmag.diff1], fvec, 8);
numFilters = 2;
filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
load([loadDir, filesep, netName], 'net');

inputData = zeros(8, numEdges, numReceivers);
[tfcomplexNN, tfcomplexNNDiff] = deal(zeros(nfft / 2, numReceivers));
input = [1; zeros(11, 1)];
windowLength = 12;
for i = 1:numReceivers
    validPath = squeeze(firstorderpathdata.diffpaths(i,:,:));
    zA = CalculateApex(rS, rR(:,i), zS, zR(:,i), wedgeLength, true);
    zA = zA(:,3);
    apexOnEdge = zA > 0 & zA < wedgeLength;
    bA = abs(thetaR(:,i) - thetaS);
    inShadow = abs(thetaR(:,i) - thetaS) > 180;
    idx = validPath & apexOnEdge & inShadow;
    pathLength = sqrt((rR(:,i) + rS) .^ 2 + (zR(:,i) - zS) .^ 2);

    rSGlobal = 2.5;
    thetaSGlobal = -45;
    source = [rSGlobal * sind(thetaSGlobal), rSGlobal * cosd(thetaSGlobal), height / 2 * ones(size(thetaSGlobal))];

    rRGlobal = 0.8;
    receiver = [rRGlobal .* sind(thetaRGlobal(i)), rRGlobal .* cosd(thetaRGlobal(i)), height / 2];
    pathLengthDir = PathLength(receiver, source);

    numInputs = sum(idx);
    if numInputs > 1

        pathLength = pathLength(idx);
        pathLengthDir = pathLengthDir .* ones(size(pathLength));
        output = CreateNNOutput(net, wedgeIndex(idx), wedgeLength(idx), thetaR(idx,i), thetaS(idx), rS(idx), rR(idx,i), zS(idx), zR(idx,i), false);

        biquad = false;
        [~, b, a] = CalculateNN(output, tfcomplex.direct(:,i), true(numInputs, 1), pathLength, nfft, fs, biquad);

        d = find(idx);
        tfcomplexNNDiffStore = zeros(nfft / 2, sum(idx));
        irLength = zeros(1, numInputs);
        irNN = cell(1, numInputs);
        for j = 1:numInputs
            %[tfmagRef, fvecRef, tfcomplexRef] = SingleWedgeInterpolated(wedgeIndex(d(j)), wedgeLength(d(j)), thetaR(d(j),i), thetaS(d(j)), rS(d(j)), rR(d(j),i), zS(d(j)), zR(d(j),i), controlparameters, false);

            [~, irNN{j}] = DelayLineIIRFilter(input, [pathLength(j), pathLengthDir(j)], windowLength, 1, b(:,:,j), a(:,:,j), c, fs, true);
            [~, tfcomplexNNDiffStore(:,j)] = IrToTf(irNN{j}, nfft);
            irLength(j) = length(irNN{j});
        end

        irStoreNN = zeros(max(irLength), numInputs);
        for j = 1:numInputs
            irStoreNN(1:irLength(j), j) = irNN{j};
        end

        irTotalNN = sum(irStoreNN, 2);
        tfcomplexNNDiff(:,i) = sum(tfcomplexNNDiffStore, 2);
        [~, diffIdx] = max(abs(tfcomplexNNDiffStore(1,:)));
        tfcomplexNNDiff(:,i) = tfcomplexNNDiffStore(:,diffIdx);
        [~, tfcomplexNNDiff(:,i)] = IrToTf(irTotalNN, nfft);

        tfcomplexNN(:,i) = tfcomplex.geom(:,i) + tfcomplexNNDiff(:,i);

%         inputData = NNInputFromGeometry(wedgeIndex(i)', wedgeLength(i)', thetaR(i,idx)', thetaS(i)', rS(i)', rR(i,idx)', zS(i)', zR(i,idx)');
%         inputData(4,:,:) = inputData(4,:,:) / 2;
%         X = dlarray(single(inputData), "CB");
% 
%         lossFunc = @()NNFilterLoss(net, X, tfmagNBand, filterFunc, false);
%         [lossNN1, ~, ~, ~, tfmagNNDiff] = dlfeval(lossFunc);
%         tfmagNN(:,i) = tfcomplex.direct(:,i) + tfcomplex.geom(:,i) + tfmagNNDiff;
    else
        [b, a] = deal(dlarray(zeros(2, 2)));
        [~, irNN] = DelayLineIIRFilter(input, [pathLength(1), pathLengthDir], windowLength, 1, b, a, c, fs, false);
        [~, tfcomplexNNDiffStore] = IrToTf(irNN, nfft);
        tfcomplexNNDiff(:,i) = tfcomplexNNDiffStore;
        tfcomplexNN(:,i) = tfcomplex.geom(:,i) + tfcomplexNNDiff(:,i);
    end
end

PlotSpectrogram(tfcomplexNN, fvec, thetaRGlobal', [-15 0], 'NN Panel', false, false, 'Receiver azimuth')
PlotSpectrogram(tfcomplexNNDiff, fvec, thetaRGlobal', [-15 0], 'NN Panel Diffraction', false, false, 'Receiver azimuth')
