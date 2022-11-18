close all

%% Input data

fs = 48e3;
nfft = 4096;
c = 344;
batchSize = 200;
controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'difforder', 1, 'saveFiles', 3);
createPlot = false;

updateRate = 200;

%% Geometry data

wedgeLength = 10;
wedgeIndex = 270;

receiverPath = [-4 2 1.6
    8 2 1.6];
source = [1 -2 1.6];
vSources = [-1 -2 1.6
    1 -2 -1.6];

[receivers, receiverDirection, receiverData, speed] = CreatePath(receiverPath, updateRate, 'scene_1', wedgeIndex, source);
numReceivers = length(receivers);

%% Create audio

windowLength = 2 * fs / updateRate;
audioLength = (numReceivers + 1) / updateRate;
[audio, audioFs] = LoopAudio('audio/whiteNoise.wav', audioLength);
audio = resample(audio,fs,audioFs)';

%% Calculate radius, theta, z

[rS, thetaS, zS] = CalculateGeometryComponents(source, wedgeIndex);
[vRS, vThetaS, vZS] = CalculateGeometryComponents(vSources, wedgeIndex);
[rR, thetaR, zR] = CalculateGeometryComponents(receivers, wedgeIndex);

%% Find apex points

% True BTM paths
[zA.ed, ~] = CalculateApex(rS, rR, zS, zR, wedgeLength);
[zA.sped, ~] = CalculateApex(vRS(2), rR, vZS(2), zR, wedgeLength);
[zA.edsp, ~] = CalculateApex(vRS(2), rR, vZS(2) + wedgeLength, zR + wedgeLength, wedgeLength);
zA.edsp = zA.edsp - wedgeLength;
[zA.spedsp, ~] = CalculateApex(rS, rR, zS + wedgeLength, zR + wedgeLength, wedgeLength);
zA.spedsp = zA.spedsp - wedgeLength;

% NN and UTD paths
[zA.NN, phii] = CalculateApex(rS, rR, zS + wedgeLength,zR + wedgeLength, 2 * wedgeLength);
zA.NN = zA.NN - wedgeLength;
[zA.vNN, vPhii] = CalculateApex(vRS(2), rR, vZS(2) + wedgeLength,zR + wedgeLength, 2 * wedgeLength);
zA.vNN = zA.vNN - wedgeLength;

%% Calculate azimuths and elevations for each path

% Diffraction paths
[azimuth, elevation] = CalculateAzimuthElevation(receiverDirection, receivers, zA);

% Direct path
[azimuth.dir, elevation.dir] = CalculateAzimuthElevation(receiverDirection, receivers, source);

% Specular paths
[azimuth.wallRef, elevation.wallRef] = CalculateAzimuthElevation(receiverDirection, receivers, vSources(1,:));
[azimuth.floorRef, elevation.floorRef] = CalculateAzimuthElevation(receiverDirection, receivers, vSources(2,:));

%% Calculate path lengths

% Direct path
pathLength.dir = PathLength(receivers, source);

% Specular paths
pathLength.wallRef = PathLength(receivers, vSources(1,:));
pathLength.floorRef = PathLength(receivers, vSources(2,:));

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

% Diffraction paths
validPath.NN = thetaR > 180 + thetaS;
validPath.vNN = thetaR > 180 + thetaS;
[validPath.ed, validPath.sped, validPath.edsp, validPath.spedsp] = deal(thetaR > 0);

%% Create delay lines

[delayedAudio, ir] = DelayLine(audio, pathLength, windowLength, validPath, c, fs);

%% BTM data

[ir.ed, ~, ~, ~, tfBtm.ed] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, createPlot);
[ir.sped, ~, ~, ~, tfBtm.sped] = SingleWedge(wedgeLength, wedgeIndex, vThetaS(2), thetaR, vRS(2), rR, vZS(2), zR, controlparameters, createPlot);
[ir.edsp, ~, ~, ~, tfBtm.edsp] = SingleWedge(wedgeLength, wedgeIndex, vThetaS(2), thetaR, vRS(2), rR, vZS(2), zR + wedgeLength, controlparameters, createPlot);
[ir.spedsp, ~, ~, ~, tfBtm.spedsp] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS + wedgeLength, zR + wedgeLength, controlparameters, createPlot);

ir.dir = ir.ed.direct;
ir.wallRef = ir.ed.geom;
ir.floorRef = ir.sped.direct;
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

% Not exact as delay then LR filter - not working need to relook.
% Previously only used for static cases.
[audioPath.utd] = DelayLineLRFilter(audio, tfmag, pathLength.NN, windowLength, validPath.NN, c, fs, nfft);
[audioPath.vUtd] = DelayLineLRFilter(audio, vTfmag, pathLength.vNN, windowLength, validPath.vNN, c, fs, nfft);

%% Continue
audioBtmTf = ConvolveTfcomplex(audio, tfBtm.ed.diff1, numReceivers);

audioBtmTest = ConvolveIR(delayedAudio.ed, irBtm, windowLength, validPath.ed);

%% HRTF

irHRTF = CreateHRTF(azimuth, elevation);
audioOut.dir = ConvolveStereoIR(delayedAudio.dir, irHRTF.dir, windowLength);
audioOut.wallRef = ConvolveStereoIR(delayedAudio.wallRef, irHRTF.wallRef, windowLength);
audioOut.floorRef = ConvolveStereoIR(delayedAudio.floorRef, irHRTF.floorRef, windowLength);
audioOut.ed = ConvolveStereoIR(audioPath.ed, irHRTF.ed, windowLength);
audioOut.sped = ConvolveStereoIR(audioPath.sped, irHRTF.sped, windowLength);
audioOut.edsp = ConvolveStereoIR(audioPath.edsp, irHRTF.edsp, windowLength);
audioOut.spedsp = ConvolveStereoIR(audioPath.spedsp, irHRTF.spedsp, windowLength);
audioOut.utd = ConvolveStereoIR(audioPath.utd, irHRTF.NN, windowLength);
audioOut.vUtd = ConvolveStereoIR(audioPath.vUtd, irHRTF.vNN, windowLength);

%% Write audio

direct = [audioOut.dir.L audioOut.dir.R];
specular = [audioOut.wallRef.L + audioOut.floorRef.L audioOut.wallRef.R + audioOut.floorRef.R];
diffracted = [audioOut.ed.L audioOut.ed.R];
specularDiffracted = [audioOut.sped.L + audioOut.edsp.L + audioOut.spedsp.L audioOut.sped.R + audioOut.edsp.R + audioOut.spedsp.R];
utd = [audioOut.utd.L + audioOut.vUtd.L audioOut.utd.R + audioOut.vUtd.R];

audiowrite('audio\bRBtm.wav', direct + specular + diffracted + specularDiffracted, fs);
audiowrite('audio\bRBtmDir.wav', direct, fs);
audiowrite('audio\bRBtmSpec.wav', specular, fs);
audiowrite('audio\bRBtmDiff.wav', diffracted, fs);
audiowrite('audio\bRBtmSpecDiff.wav', specularDiffracted, fs);
audiowrite('audio\bRUtd.wav', direct + specular + utd, fs);
audiowrite('audio\bRUtdDiff.wav', utd, fs);

direct = audioPath.dir;
specular = audioPath.wallRef + audioPath.floorRef;
diffracted = audioPath.ed;
specularDiffracted = audioPath.sped + audioPath.edsp + audioPath.spedsp;
utd = audioPath.utd + audioPath.vUtd;

audiowrite('audio\Btm.wav', direct + specular + diffracted + specularDiffracted, fs);
audiowrite('audio\BtmDir.wav', direct, fs);
audiowrite('audio\BtmSpec.wav', specular, fs);
audiowrite('audio\BtmDiff.wav', diffracted, fs);
audiowrite('audio\BtmSpecDiff.wav', specularDiffracted, fs);
audiowrite('audio\Utd.wav', direct + specular + utd, fs);
audiowrite('audio\UtdDiff.wav', utd, fs);

audiowrite('audio\trueBTMtfwn.wav', audioBtmTf, fs)
audiowrite('audio\delayedBTMwn.wav', audioBtmTest, fs)

%% Figures

close all

PlotSpectrogramOfWAV('audio\bRBtm.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\bRBtmDir.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\bRBtmSpec.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\bRBtmDiff.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\bRBtmSpecDiff.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\bRUtd.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\bRUtdDiff.wav', [-70 0], nfft);

PlotSpectrogramOfWAV('audio\Btm.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\BtmDir.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\BtmSpec.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\BtmDiff.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\BtmSpecDiff.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\Utd.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\UtdDiff.wav', [-70 0], nfft);

PlotSpectrogramOfWAV('audio\trueBTMwn.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\trueBTMtfwn.wav', [-70 0], nfft);
PlotSpectrogramOfWAV('audio\delayedBTMwn.wav', [-70 0], nfft);

input = [1; zeros(5, 1)];

createPlot = false;
arg = cosd(90 - phii);
l = rS / arg;
m = rR / arg;

time = numReceivers / updateRate;

[audio, fs] = audioread('audio\whiteNoise.wav');

ir = CreateHRTF(azimuth, elevation);
audioDir = ConvolveStereoIR(audio, ir.Dir, windowLength);
audioGeom = ConvolveStereoIR(audio, ir.Geom, windowLength);



load 'ReferenceHRTF.mat' hrtfData sourcePosition

hrtfData = permute(double(hrtfData),[2,3,1]);

sourcePosition = sourcePosition(:,[1,2]);

n = 100;
durationPerPosition = 0.2;
%durationPerPosition = 1;
desiredAz = [-120;-60;0;60;120;0;-120;120];
desiredAz = linspace(90,90,n)';
desiredEl = [-90;90;45;0;-45;0;45;45];
desiredEl = linspace(0,180,n)';

desiredPosition = [desiredAz desiredEl];

interpolatedIR  = interpolateHRTF(hrtfData,sourcePosition,desiredPosition);

leftIR = squeeze(interpolatedIR(:,1,:));
rightIR = squeeze(interpolatedIR(:,2,:));

desiredFs = 48e3;
audioFilePath = 'audio\imagination.wav';
[audio, fs] = LoopAudio(audioFilePath, 0.8 * n * durationPerPosition);
audio = 0.8*resample(audio,desiredFs,fs);
audioFilePath = 'audio\imagination48k.wav';
audiowrite(audioFilePath,audio,desiredFs);

audioLength = length(audio);

fileReader = dsp.AudioFileReader(audioFilePath);
deviceWriter = audioDeviceWriter('SampleRate',fileReader.SampleRate);

leftFilter = dsp.FIRFilter('NumeratorSource','Input port');
rightFilter = dsp.FIRFilter('NumeratorSource','Input port');

samplesPerPosition = durationPerPosition*fileReader.SampleRate;
samplesPerPosition = samplesPerPosition - rem(samplesPerPosition,fileReader.SamplesPerFrame);

sourcePositionIndex = 1;
samplesRead = 0;
while ~isDone(fileReader)
    audioIn = fileReader();
    samplesRead = samplesRead + fileReader.SamplesPerFrame;
    
    leftChannel = leftFilter(audioIn,leftIR(sourcePositionIndex,:));
    rightChannel = rightFilter(audioIn,rightIR(sourcePositionIndex,:));
    
    deviceWriter([leftChannel,rightChannel]);
    
    samplesRemaining = audioLength - samplesRead;
    if sourcePositionIndex > 98
        x = 2;
    end
    if samplesRemaining < 1000
        x = 1;
    end
    if mod(samplesRead,samplesPerPosition) == 0
        sourcePositionIndex = sourcePositionIndex + 1;
    end
end

release(deviceWriter)
release(fileReader)