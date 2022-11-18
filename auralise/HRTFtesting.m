clear all
close all

%% Input data

fs = 48e3;
nfft = 4096;
c = 344;
batchSize = 200;
controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'difforder', 1, 'saveFiles', true);

updateRate = 100;

%% Geometry data

wedgeLength = 10;
wedgeIndex = 270;

receiverPath = [-4 2 1.6
    1 2 1.6
    8 -4 1.6];
source = [1 -2 1.6];
vSources = [-1 -2 1.6
    1 -2 -1.6];

[receivers, receiverDirection, receiverData, speed] = CreatePath(receiverPath, updateRate, 'scene_1', wedgeIndex, source);
numReceivers = length(receivers);

%% Create audio

windowLength = 2 * fs / updateRate;
audioLength = numReceivers / updateRate;
[audio, audioFs] = LoopAudio('audio/whiteNoise.wav', audioLength);
audio = resample(audio,fs,audioFs);

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
pathLength.vNN = DiffractionPathLength(receivers, vSources(2,:), zA.VNN);

%% Calculate valid paths

% Direct path
validPath.dir = thetaR <= 180;

% Specular paths
validPath.wallRef = thetaR <= 180 - 2 * thetaS;
validPath.floorRef = thetaR <= 180;

% Diffraction paths
validPath.NN = thetaR > 180;
validPath.vNN = thetaR > 180;
[validPath.ed, validPath.sped, validPath.edsp, validPath.spedsp] = deal(thetaR > 0);

%% Create delay lines

[delayedAudio, ir] = DelayLine(audio, pathLength, windowLength, validPath, c, fs);
DelayLineBiquad

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