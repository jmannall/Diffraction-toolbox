close all
clear all

%% Input data

fs = 48e3;
nfft = 4096;
c = 344;
batchSize = 200;
controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'difforder', 1, 'saveFiles', 3);
createPlot = false;

speed = 1.3; % Average walking speed

% Minimum update rate so that the delay length changes by < 1 sample each
% update

F = factor(fs);

allFactors = [];
for i = 1:length(F)
    allFactors = [allFactors; unique(prod(nchoosek(F, i), 2))];
end
allFactors = sort(unique(allFactors));

frameRate = 30;
% Find lowest update rate that is often enough to prevent delay shifts
% greater than 1 sample and a multiple of fs.
updateRate = ceil(speed * fs / c);
updateRate = allFactors(find(updateRate < allFactors, 1));

%% Geometry data

wedgeLength = 10;
wedgeIndex = 270;

receiverPath = [-4 2 1.6
    8 2 1.6];
receiverHeading = [0 -1 0];
source = [1 -2 1.6];
vSources = [-1 -2 1.6
    1 -2 -1.6
    -1 -2 -1.6];

scene = 'scene_1';
[receivers, receiverDirection, receiverHeading, receiverData, speed] = CreatePath(receiverPath, updateRate, frameRate, scene, wedgeIndex, source, receiverHeading, speed);
numReceivers = length(receivers);

%% Create audio

windowLength = 2 * fs / updateRate;
audioLength = (numReceivers + 1) / updateRate;
[audio, audioFs] = LoopAudio('audio/whiteNoise.wav', audioLength);
audio = resample(audio,fs,audioFs)';

%% LR testing

% pathLength = 0.001 * ones(1, numReceivers);
% validPath = true(size(pathLength));
% tfmag = [0 -6 -18 -12];
% 
% audiowrite(['audio', filesep, 'orginal.wav'], audio, fs);
% PlotSpectrogramOfWAV(['audio', filesep, 'orginal.wav'], [-70 0], nfft);
% 
% %% Delay line LR
% output = DelayLineLR(audio, pathLength, windowLength, validPath, tfmag, c, fs);
% 
% audiowrite(['audio', filesep, 'LRfilter.wav'], output, fs);
% PlotSpectrogramOfWAV(['audio', filesep, 'LRfilter.wav'], [-70 0], nfft);

%% Calculate radius, theta, z

[rS, thetaS, zS] = CalculateGeometryComponents(source, wedgeIndex);
[vRS, vThetaS, vZS] = CalculateGeometryComponents(vSources, wedgeIndex);
[rR, thetaR, zR] = CalculateGeometryComponents(receivers, wedgeIndex);

%% Find apex points

% True BTM paths
[zA.ed, ~] = CalculateApex(rS, rR, zS, zR, wedgeLength);
[zA.sped, ~] = CalculateApex(vRS(2), rR, vZS(2), zR, wedgeLength);
[zA.edsp, ~] = CalculateApex(vRS(2), rR, vZS(2) + wedgeLength, zR + wedgeLength, wedgeLength);
zA.edsp(:,3) = wedgeLength - zA.edsp(:,3);
[zA.spedsp, ~] = CalculateApex(rS, rR, zS + wedgeLength, zR + wedgeLength, wedgeLength);
zA.spedsp(:,3) = zA.spedsp(:,3) - wedgeLength;

% NN and UTD paths
[zA.NN, phii] = CalculateApex(rS, rR, zS + wedgeLength,zR + wedgeLength, 2 * wedgeLength);
zA.NN(:,3) = zA.NN(:,3) - wedgeLength;
[zA.vNN, vPhii] = CalculateApex(vRS(2), rR, vZS(2) + wedgeLength,zR + wedgeLength, 2 * wedgeLength);
zA.vNN(:,3) = zA.vNN(:,3) - wedgeLength;

%% Calculate azimuths and elevations for each path

% Diffraction paths
[azimuth, elevation] = CalculateAzimuthElevation(receiverHeading, receivers, zA);

% Direct path
[azimuth.dir, elevation.dir] = CalculateAzimuthElevation(receiverHeading, receivers, source);

% Specular paths
[azimuth.wallRef, elevation.wallRef] = CalculateAzimuthElevation(receiverHeading, receivers, vSources(1,:));
[azimuth.floorRef, elevation.floorRef] = CalculateAzimuthElevation(receiverHeading, receivers, vSources(2,:));
[azimuth.floorWallRef, elevation.floorWallRef] = CalculateAzimuthElevation(receiverHeading, receivers, vSources(3,:));

x = thetaR - thetaS;
figure
plot(x, azimuth.dir)
hold on
plot(x, azimuth.wallRef)
plot(x, azimuth.ed)
legend('direct', 'specular', 'diffraction')
title('azimuth')

% figure
% plot(x, elevation.dir)
% hold on
% plot(x, elevation.wallRef)
% plot(x, elevation.ed)
% legend('direct', 'specular', 'diffraction')
% title('elevation')

%% Calculate path lengths

% Direct path
pathLength.dir = PathLength(receivers, source);

% Specular paths
pathLength.wallRef = PathLength(receivers, vSources(1,:));
pathLength.floorRef = PathLength(receivers, vSources(2,:));
pathLength.floorWallRef = PathLength(receivers, vSources(3,:));

% Diffraction paths
pathLength.ed = DiffractionPathLength(receivers, source, zA.ed);
pathLength.sped = DiffractionPathLength(receivers, vSources(2,:), zA.sped);
pathLength.edsp = DiffractionPathLength(receivers, vSources(2,:), zA.edsp);
pathLength.spedsp = DiffractionPathLength(receivers, source, zA.spedsp);

pathLength.NN = DiffractionPathLength(receivers, source, zA.NN);
pathLength.vNN = DiffractionPathLength(receivers, vSources(2,:), zA.vNN);

%% Calculate valid paths

% Direct path
validPath.dir = thetaR <= 180 + thetaS;

% Specular paths
validPath.wallRef = thetaR <= 180 - thetaS;
validPath.floorRef = thetaR <= 180 + thetaS;
validPath.floorWallRef = thetaR <= 180 - thetaS;

% Diffraction paths
validPath.NN = thetaR > 180 + thetaS;
validPath.vNN = thetaR > 180 + thetaS;
[validPath.ed, validPath.sped, validPath.edsp, validPath.spedsp] = deal(thetaR > 0);

%% Create delay lines

[delayedAudio, ir] = DelayLine(audio, pathLength, windowLength, validPath, c, fs);

%% BTM data

[ir.ed, ~, ~, ~, tfBtm.ed] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, createPlot);
[ir.sped, ~, ~, ~, tfBtm.sped] = SingleWedge(wedgeLength, wedgeIndex, vThetaS(2), thetaR, vRS(2), rR, vZS(2), zR, controlparameters, createPlot);
[ir.edsp, ~, ~, ~, tfBtm.edsp] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR - 2 * zR, controlparameters, createPlot);
[ir.spedsp, ~, ~, fvec, tfBtm.spedsp] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS + wedgeLength, zR + wedgeLength, controlparameters, createPlot);

% ir.dir = ir.ed.direct;
% ir.wallRef = ir.ed.geom;
% ir.floorRef = ir.sped.direct;
ir.ed = ir.ed.diff1;
ir.sped = ir.sped.diff1;
ir.edsp = ir.edsp.diff1;
ir.spedsp = ir.spedsp.diff1;

[idxStart, idxEnd, irLength] = deal(zeros(numReceivers, 1));
for i = 1:numReceivers
    idxStart(i) = find(ir.ed(:,i), 1, "first");
    idxEnd(i) = find(ir.ed(:,i), 1, "last");
    irLength(i) = idxEnd(i) - idxStart(i) + 1;
end
nDir = max(irLength);
irBtm = zeros(nDir, numReceivers);
for i = 1:numReceivers
    idx = idxStart(i):idxEnd(i);
    irBtm(1:irLength(i),i) = pathLength.ed(i) * ir.ed(idx,i);
end

%% BTM audio

% audioPath.dir = ConvolveIR(audio, ir.dir, windowLength, validPath.dir);
% audioPath.wallRef = ConvolveIR(audio, ir.wallRef, windowLength, validPath.wallRef);
% audioPath.floorRef = ConvolveIR(audio, ir.floorRef, windowLength, validPath.floorRef);
audioPath.ed = ConvolveIR(audio, ir.ed, windowLength, validPath.ed);
audioPath.sped = ConvolveIR(audio, ir.sped, windowLength, validPath.sped);
audioPath.edsp = ConvolveIR(audio, ir.edsp, windowLength, validPath.edsp);
audioPath.spedsp = ConvolveIR(audio, ir.spedsp, windowLength, validPath.spedsp);

%% UTD audio

[tfmag, vTfmag] = deal(zeros(numReceivers, 4));
for i = 1:numReceivers
    [tfmag(i,:), ~, ~] = SingleUTDWedgeInterpolated(thetaS, thetaR(i), rS, rR(i), wedgeIndex, phii(i), controlparameters);
    [vTfmag(i,:), ~, ~] = SingleUTDWedgeInterpolated(vThetaS(2), thetaR(i), vRS(2), rR(i), wedgeIndex, vPhii(i), controlparameters);
end

[audioPath.utd, ir.utd] = DelayLineLR(audio, pathLength.NN, windowLength, validPath.ed, tfmag, c, fs);
[audioPath.vUtd, ir.vUtd] = DelayLineLR(audio, pathLength.vNN, windowLength, validPath.ed, vTfmag, c, fs);

%% NN audio

% loadPath = ['NNSaves', filesep, 'Biquad-700_00060474-14515-089849-080573-3'];
loadPath = ['NNSaves', filesep, 'IIR-20000_0001-1-09-099-2.mat'];
load(loadPath, "net");

% X = [wedgeLength, wedgeIndex, thetaR - thetaS]  ;  
% Where one corresponds to the shortest radius. And the z axis zero is the
% edge corner nearest to zOne
bendingAngle = thetaR - thetaS;
minAngle = min(thetaS, wedgeIndex - thetaR);


sourceIsROne = rS < rR;
rOne(sourceIsROne, 1) = rS;
rOne(~sourceIsROne, 1) = rR(~sourceIsROne);
rTwo(sourceIsROne, 1) = rR(sourceIsROne);
rTwo(~sourceIsROne, 1) = rS;

zOne(sourceIsROne, 1) = zS + wedgeLength;
zOne(~sourceIsROne, 1) = zR(~sourceIsROne, 1) + wedgeLength;
zTwo(sourceIsROne, 1) = zR(sourceIsROne, 1) + wedgeLength;
zTwo(~sourceIsROne, 1) = zS + wedgeLength;

longWedgeLength = wedgeLength * 2;

flipCorner = zOne > wedgeLength;
zOne(flipCorner, 1) = wedgeLength - zOne(flipCorner, 1);
zTwo(flipCorner, 1) = wedgeLength - zTwo(flipCorner, 1);

const = ones(numReceivers, 1);
X = ([deg2rad(wedgeIndex) * const, deg2rad(bendingAngle), deg2rad(minAngle), longWedgeLength * const, rOne, rTwo, zOne, zTwo])';
X = dlarray(single(X), "CB");
output = predict(net, X);

tfmag = zeros(nfft / 2, numReceivers);
tfcomplexNNRef = zeros(nfft / 2, numReceivers);
for i = 1:numReceivers
    [tfmagTest(:,i), ~, tfcomplexNNRef(:,i)] = SingleWedgeInterpolated(longWedgeLength,wedgeIndex,thetaS,thetaR(i),rS,rR(i),zS,zR(i),controlparameters,false);
end

tfcomplexNNRef(:, validPath.dir) = tfcomplexNNRef(:, validPath.dir) ./ pathLength.dir(validPath.dir)';
tfcomplexNNRef(:, validPath.NN) = tfcomplexNNRef(:, validPath.NN) ./ pathLength.NN(validPath.NN)';

PlotSpectrogram(tfcomplexNNRef, fvec, thetaR - thetaS, [-70 0], 'NN', false, false, 'Bending Angle')

%% go
numFilters = 2;
biquad = false;
[targetData, fc, fidx] = CreateFrequencyNBands(tfmagTest, fvec, 12);

if biquad
    loss = BiquadLoss(output, targetData, numFilters, nfft, fs, fidx);

    [zReal, zImg, pReal, pImg, k] = CreateZPKFromNNOutput(output, numFilters);

    [b, a] = BiquadCoefficients(zReal, zImg, pReal, pImg, k, numFilters, numReceivers);
    [tfmagNN, ~, tfcomplexNNstore] = CreateBiquad(zReal, zImg, pReal, pImg, k, nfft, fs);
else
    loss = IIRFilterLoss(output, targetData, numFilters, nfft, fs, fidx);

    [z, p, k] = CreateIIRFromNNOutput(output, numFilters);
    [b, a] = IIRFilterCoefficients(z, p, k, numFilters, numReceivers);
    [tfmagNN, ~, tfcomplexNNstore] = CreateIIRFilter(z, p, k, nfft, fs);
end

[tfmag, tfcomplex] = IrToTf(ir, nfft);

tfcomplexNN = tfcomplex.dir;
tfcomplexNN(:,validPath.NN) = extractdata(tfcomplexNNstore(:,validPath.NN) ./ pathLength.NN(validPath.NN)');
PlotSpectrogram(tfcomplexNN, fvec, thetaR - thetaS, [-70 0], 'NN', false, false, 'Bending Angle')

idx = 50;
figure
semilogx(fvec, tfmagTest(:,idx))
hold on
semilogx(fvec, extractdata(tfmagNN(:,idx)))
xlim([20 20e3])

[audioPath.NN, ir.NN] = DelayLineBiquad(audio, pathLength.NN, windowLength, validPath.ed, b, a, numBiquads, c, fs, validPath.NN);

%% Transfer functions

[tfmag, tfcomplex] = IrToTf(ir, nfft);

close all

PlotSpectrogram([tfcomplex.dir], fvec, thetaR - thetaS, [-70 0], 'Direct', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.wallRef], fvec, thetaR - thetaS, [-70 0], 'Wall reflection', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.floorRef], fvec, thetaR - thetaS, [-70 0], 'Floor Reflection', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.floorWallRef], fvec, thetaR - thetaS, [-70 0], 'Floor Wall Reflection', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.ed], fvec, thetaR - thetaS, [-70 0], 'Edge Diffraction', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.sped], fvec, thetaR - thetaS, [-70 0], 'Floor Edge Diffraction', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.edsp], fvec, thetaR - thetaS, [-70 0], 'Edge Floor Diffarction', false, false, 'Bending Angle')
PlotSpectrogram([tfcomplex.spedsp], fvec, thetaR - thetaS, [-70 0], 'Floor Edge Floor Diffraction', false, false, 'Bending Angle')

geomField = [tfcomplex.dir] + [tfcomplex.wallRef] + [tfcomplex.floorRef] + [tfcomplex.floorWallRef];
diffFieldBtm = [tfcomplex.ed] + [tfcomplex.sped] + [tfcomplex.edsp] + [tfcomplex.spedsp];
completeFieldBtm = geomField + diffFieldBtm;
PlotSpectrogram(geomField, fvec, thetaR - thetaS, [-70 0], 'Geometric Field', false, false, 'Bending Angle')
PlotSpectrogram(diffFieldBtm, fvec, thetaR - thetaS, [-70 0], 'BTM Diffracted Field', false, false, 'Bending Angle')
PlotSpectrogram(completeFieldBtm, fvec, thetaR - thetaS, [-70 0], 'BTM Complete Field', false, false, 'Bending Angle')

diffFieldUtd = [tfcomplex.utd] + [tfcomplex.vUtd];
completeFieldUtd = [tfcomplex.wallRef] + [tfcomplex.floorWallRef] + diffFieldUtd;
PlotSpectrogram(diffFieldUtd, fvec, thetaR - thetaS, [-70 0], 'UTD Diffracted Field', false, false, 'Bending Angle')
PlotSpectrogram(completeFieldUtd, fvec, thetaR - thetaS, [-70 0], 'UTD Complete Field', false, false, 'Bending Angle')

PlotSpectrogram([tfcomplex.NN], fvec, thetaR - thetaS, [-70 0], 'Single NN Path', false, false, 'Bending Angle')

%% HRTF

irHRTF = CreateHRTF(azimuth, elevation);
audioOut.dir = ConvolveStereoIR(delayedAudio.dir, irHRTF.dir, windowLength);
audioOut.wallRef = ConvolveStereoIR(delayedAudio.wallRef, irHRTF.wallRef, windowLength);
audioOut.floorRef = ConvolveStereoIR(delayedAudio.floorRef, irHRTF.floorRef, windowLength);
audioOut.floorWallRef = ConvolveStereoIR(delayedAudio.floorWallRef, irHRTF.floorWallRef, windowLength);
audioOut.ed = ConvolveStereoIR(audioPath.ed, irHRTF.ed, windowLength);
audioOut.sped = ConvolveStereoIR(audioPath.sped, irHRTF.sped, windowLength);
audioOut.edsp = ConvolveStereoIR(audioPath.edsp, irHRTF.edsp, windowLength);
audioOut.spedsp = ConvolveStereoIR(audioPath.spedsp, irHRTF.spedsp, windowLength);
audioOut.utd = ConvolveStereoIR(audioPath.utd, irHRTF.NN, windowLength);
audioOut.vUtd = ConvolveStereoIR(audioPath.vUtd, irHRTF.vNN, windowLength);

%% Write audio

geometric = SumAudioOut(audioOut, {'dir', 'wallRef', 'floorRef', 'floorWallRef'});
geometricUtdNN = SumAudioOut(audioOut, {'wallRef', 'floorWallRef'});
diffractedBtm = SumAudioOut(audioOut, {'ed', 'sped', 'edsp', 'spedsp'});
diffractedUtd = SumAudioOut(audioOut, {'utd', 'vUtd'});

audioFilePath = ['audio\', scene];

audiowrite([audioFilePath, '_bRBtm.wav'], geometric + diffractedBtm, fs);
audiowrite([audioFilePath, '_bRGeometric.wav'], geometric, fs);
audiowrite([audioFilePath, '_bRUtd.wav'], geometricUtdNN + diffractedUtd, fs);
audiowrite([audioFilePath, '_bRUtdNoWallRef.wav'], diffractedUtd, fs);

geometric = delayedAudio.dir + delayedAudio.wallRef + delayedAudio.floorRef + delayedAudio.floorWallRef;
geometricUtdNN = delayedAudio.wallRef + delayedAudio.floorWallRef;
diffractedBtm = audioPath.ed + audioPath.sped + audioPath.edsp + audioPath.spedsp;
diffractedUtd = audioPath.utd + audioPath.vUtd;

audiowrite([audioFilePath, '_Btm.wav'], geometric + diffractedBtm, fs);
audiowrite([audioFilePath, '_Geometric.wav'], geometric, fs);
audiowrite([audioFilePath, '_Utd.wav'], geometricUtdNN + diffractedUtd, fs);
audiowrite([audioFilePath, '_UtdNoWallRef.wav'], diffractedUtd, fs);

%% Figures

close all

PlotSpectrogramOfWAV([audioFilePath, '_bRBtm.wav'], [-70 0], nfft);
PlotSpectrogramOfWAV([audioFilePath, '_bRGeometric.wav'], [-70 0], nfft);
PlotSpectrogramOfWAV([audioFilePath, '_bRUtd.wav'], [-70 0], nfft);
PlotSpectrogramOfWAV([audioFilePath, '_bRUtdNoWallRef.wav'], [-70 0], nfft);

PlotSpectrogramOfWAV([audioFilePath, '_Btm.wav'], [-70 0], nfft);
PlotSpectrogramOfWAV([audioFilePath, '_Geometric.wav'], [-70 0], nfft);
PlotSpectrogramOfWAV([audioFilePath, '_Utd.wav'], [-70 0], nfft);
PlotSpectrogramOfWAV([audioFilePath, '_UtdNoWallRef.wav'], [-70 0], nfft);