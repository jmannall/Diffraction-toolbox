close all
clear all

wallThickness = 1;
wallHeight = 10;
bendingAngle = 5:360;
minAngle = 30;
[radiusR, radiusS] = deal(1);
[zS, zR] = deal(wallHeight / 2);

nfft = 4096;
fs = 48000;
controlparameters = struct('nfft', nfft, 'fs', fs, 'difforder', 3);

geometry = GeometryWall(wallThickness, bendingAngle, minAngle, radiusR);
bendingAngle = geometry.bendingAngle;
[result, geometry] = SingleWallArray(geometry, wallHeight, radiusS, radiusR, zS, zR, controlparameters);

ir = [result.ir];
tfcomplex = [result.tfcomplex];
fvec = result(1).fvec;
PlotSpectogram([tfcomplex.complete], fvec, bendingAngle, [-70 0], 'BTM wall all', false, 'Bending Angle')
PlotSpectogram([tfcomplex.direct], fvec, bendingAngle, [-70 0], 'BTM wall direct', false, 'Bending Angle')
PlotSpectogram([tfcomplex.geom], fvec, bendingAngle, [-70 0], 'BTM wall specular', false, 'Bending Angle')
PlotSpectogram([tfcomplex.diff1], fvec, bendingAngle, [-70 0], 'BTM wall diff1', true, 'Bending Angle')
PlotSpectogram([tfcomplex.diff2], fvec, bendingAngle, [-70 0], 'BTM wall diff2', true, 'Bending Angle')
PlotSpectogram([tfcomplex.diff3], fvec, bendingAngle, [-70 0], 'BTM wall diff3', true, 'Bending Angle')
PlotSpectogram([tfcomplex.direct] + [tfcomplex.geom] + [tfcomplex.diff1], fvec, bendingAngle, [-70 0], 'BTM wall up to diff1', false, 'Bending Angle')
PlotSpectogram([tfcomplex.direct] + [tfcomplex.geom] + [tfcomplex.diff1] + [tfcomplex.diff2], fvec, bendingAngle, [-70 0], 'BTM wall up to diff2', false, 'Bending Angle')


%% audio

fileStem = 'audio\whiteNoise';
audioFilePath = [fileStem, '.wav'];
audio = NormalizeAudio(audioFilePath, nfft);

exampleLength = 7.25;
[audio, fs, audioLength] = LoopAudio([fileStem, '_Normalised.wav'], exampleLength);

updateRate = 20;
samplesPerUpdate = fs / updateRate;
numWindows = 2 * floor(audioLength / samplesPerUpdate) - 2;

validPath = bendingAngle > 0;
test = [ir.complete];
BTM = ConvolveIR(audio, [ir.complete], validPath, samplesPerUpdate);

savePath = [fileStem, '_BTMWallAll.wav'];
audiowrite(savePath, BTM, fs)
PlotSpectogramOfWAV([fileStem, '_BTMWallAll.wav'], [-70 0], nfft);
